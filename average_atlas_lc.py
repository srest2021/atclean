#!/usr/bin/env python
"""
@author: Sofia Rest
"""

import sys, argparse, configparser, re
from copy import deepcopy
import pandas as pd
import numpy as np

import sigmacut
from pdastro import pdastrostatsclass, AandB, AnotB
from atlas_lc import atlas_lc
from plot_atlas_lc import plot_atlas_lc

class average_atlas_lc():
	def __init__(self):
		# input/output
		self.input_dir = None
		self.output_dir = None
		self.overwrite = True

		# flags
		self.flags = {'flag_chisquare':0x1,

					  'flag_uncertainty':0x2,

					  'flag_controls_bad':0x400000,
					  'flag_controls_questionable':0x80000,
					  'flag_controls_x2':0x100,
					  'flag_controls_stn':0x200,
					  'flag_controls_Nclip':0x400,
					  'flag_controls_Ngood':0x800,

					  'flag_badday':0x800000,
					  'flag_ixclip':0x1000,
					  'flag_smallnum':0x2000}

		self.mjd_bin_size = None
		self.keep_empty_bins = False
		self.flux2mag_sigmalimit = None

		self.Nclip_max = None
		self.Ngood_min = None 
		self.x2_max = None

	# define command line arguments
	def define_args(self, parser=None, usage=None, conflict_handler='resolve'):
		if parser is None:
			parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)
		
		parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')
		parser.add_argument('-m', '--mjd_bin_size', type=float, default=None, help='MJD bin size in days for averaging')

		parser.add_argument('-p', '--plot', default=False, action='store_true', help='plot each cut and save into PDF file')
		parser.add_argument('--xlim_lower', type=float, default=None, help='if plotting, manually set lower x axis limit to a certain MJD')
		parser.add_argument('--xlim_upper', type=float, default=None, help='if plotting, manually set upper x axis limit to a certain MJD')
		parser.add_argument('--ylim_lower', type=float, default=None, help='if plotting, manually set lower y axis limit to a certain uJy')
		parser.add_argument('--ylim_upper', type=float, default=None, help='if plotting, manually set upper y axis limit to a certain uJy')

		parser.add_argument('-f','--cfg_filename', default='atlas_lc_settings.ini', type=str, help='file name of ini file with settings for this class')
		parser.add_argument('--dont_overwrite', default=False, action='store_true', help='don\'t overwrite existing file with same file name')

		return parser
	
	# load config settings from file and reconcile with command arguments
	def load_settings(self, args):
		print('LOADING SETTINGS FROM CONFIG FILE AND CMD ARGUMENTS...')

		cfg = configparser.ConfigParser()
		try:
			print(f'Loading config file at {args.cfg_filename}')
			cfg.read(args.cfg_filename)
		except Exception as e:
			raise RuntimeError(f'ERROR: Could not load config file at {args.cfg_filename}!')

		self.input_dir = cfg['Input/output settings']['input_dir']
		print(f'Light curve .txt files input directory: {self.input_dir}')
		self.output_dir = cfg['Input/output settings']['output_dir']
		print(f'Light curve .txt files output directory: {self.output_dir}')
		self.overwrite = not args.dont_overwrite
		print(f'Overwrite existing light curve files: {self.overwrite}')
		self.flux2mag_sigmalimit = cfg['Input/output settings']['flux2mag_sigmalimit']
		print(f'Sigma limit when converting flux to magnitude (magnitudes are limits when dmagnitudes are NaN): {self.flux2mag_sigmalimit}')
		print(f'Plotting: {args.plot}')

		self.mjd_bin_size = args.mjd_bin_size if not(args.mjd_bin_size is None) else float(cfg['Averaging settings']['mjd_bin_size'])
		print(f'MJD bin size: {self.mjd_bin_size} days')
		self.keep_empty_bins = bool(cfg['Averaging settings']['keep_empty_bins'])
		print(f'Keep empty bins and store as NaN in averaged light curve: {self.keep_empty_bins}')
		
		self.Nclip_max = int(cfg['Averaging settings']['Nclip_max'])
		self.Ngood_min = int(cfg['Averaging settings']['Ngood_min'])
		self.x2_max = float(cfg['Averaging settings']['x2_max'])
		print(f'MJD bin bounds for not flagging as bad day: ')
		print(f'# Maximum number of clipped measurements (Nclip_max): {self.Nclip_max}')
		print(f'# Minimum number of good measurements (Ngood_min): {self.Ngood_min}')
		print(f'# Maximum chi-square (x2_max): {self.x2_max}')

	# drop mask column and any extra/unneeded columns from previous iterations
	def drop_extra_columns(self, lc):
		dropcols=[]

		for col in ['Noffsetlc', '__tmp_SN']:
			if col in lc.pdastro.t.columns:
				dropcols.append(col)
		for col in lc.pdastro.t.columns:
			if re.search('^c\d_',col): 
				dropcols.append(col)

		# drop any extra columns
		if len(dropcols)>0: 
			print('Dropping extra columns: ',dropcols)
			lc.pdastro.t.drop(columns=dropcols,inplace=True)

		return lc

	def average_lc(self, lc, avglc):
	    mjd = int(np.amin(lc.pdastro.t['MJD']))
	    mjd_max = int(np.amax(lc.pdastro.t['MJD']))+1

	    good_i = lc.pdastro.ix_unmasked('Mask', maskval=self.flags['flag_chisquare']|self.flags['flag_uncertainty']|self.flags['flag_controls_bad'])

	    while mjd <= mjd_max:
	        range_i = lc.pdastro.ix_inrange(colnames=['MJD'], lowlim=mjd, uplim=mjd+self.mjd_bin_size, exclude_uplim=True)
	        range_good_i = AandB(range_i,good_i)

	        # add new row to averaged light curve if keep_empty_bins or any measurements present
	        if self.keep_empty_bins or len(range_i) >= 1:
	            new_row = {'MJDbin':mjd+0.5*self.mjd_bin_size, 'Nclip':0, 'Ngood':0, 'Nexcluded':len(range_i)-len(range_good_i), 'Mask':0}
	            avglc_index = avglc.pdastro.newrow(new_row)
	        
	        # if no measurements present, flag or skip over day
	        if len(range_i) < 1:
	            if self.keep_empty_bins:
	                avglc.update_mask_col(self.flags['flag_badday'], [avglc_index])
	            mjd += self.mjd_bin_size
	            continue
	        
	        # if no good measurements, average values anyway and flag
	        if len(range_good_i) < 1:
	            # average flux
	            lc.pdastro.calcaverage_sigmacutloop('uJy', noisecol='duJy', indices=range_i, Nsigma=3.0, median_firstiteration=True)
	            fluxstatparams = deepcopy(lc.pdastro.statparams)

	            # get average mjd
	            # TO DO: SHOULD NOISECOL HERE BE DUJY OR NONE??
	            lc.pdastro.calcaverage_sigmacutloop('MJD', noisecol='duJy', indices=fluxstatparams['ix_good'], Nsigma=0, median_firstiteration=False)
	            avg_mjd = lc.pdastro.statparams['mean']

	            # add row and flag
	            avglc.pdastro.add2row(avglc_index, {'MJD':avg_mjd, 
	                                                'uJy':fluxstatparams['mean'], 
	                                                'duJy':fluxstatparams['mean_err'], 
	                                                'stdev':fluxstatparams['stdev'],
	                                                'x2':fluxstatparams['X2norm'],
	                                                'Nclip':fluxstatparams['Nclip'],
	                                                'Ngood':fluxstatparams['Ngood'],
	                                                'Mask':0})
	            lc.update_mask_col(self.flags['flag_badday'], range_i)
	            avglc.update_mask_col(self.flags['flag_badday'], [avglc_index])

	            mjd += self.mjd_bin_size
	            continue
	        
	        # average good measurements
	        lc.pdastro.calcaverage_sigmacutloop('uJy', noisecol='duJy', indices=range_good_i, Nsigma=3.0, median_firstiteration=True)
	        fluxstatparams = deepcopy(lc.pdastro.statparams)

	        if fluxstatparams['mean'] is None or len(fluxstatparams['ix_good']) < 1:
	            lc.update_mask_col(self.flags['flag_badday'], range_i)
	            avglc.update_mask_col(self.flags['flag_badday'], [avglc_index])
	            mjd += self.mjd_bin_size
	            continue

	        # get average mjd
	        # TO DO: SHOULD NOISECOL HERE BE DUJY OR NONE??
	        lc.pdastro.calcaverage_sigmacutloop('MJD', noisecol='duJy', indices=fluxstatparams['ix_good'], Nsigma=0, median_firstiteration=False)
	        avg_mjd = lc.pdastro.statparams['mean']

	        # add row to averaged light curve
	        avglc.pdastro.add2row(avglc_index, {'MJD':avg_mjd, 
	                                            'uJy':fluxstatparams['mean'], 
	                                            'duJy':fluxstatparams['mean_err'], 
	                                            'stdev':fluxstatparams['stdev'],
	                                            'x2':fluxstatparams['X2norm'],
	                                            'Nclip':fluxstatparams['Nclip'],
	                                            'Ngood':fluxstatparams['Ngood'],
	                                            'Mask':0})
	        
	        # flag clipped measurements in lc
	        if len(fluxstatparams['ix_clip']) > 0:
	            lc.update_mask_col(self.flags['flag_ixclip'], fluxstatparams['ix_clip'])
	        
	        # if small number within this bin, flag measurements
	        if len(range_good_i) < 3:
	            lc.update_mask_col(self.flags['flag_smallnum'], range_good_i) # TO DO: CHANGE TO RANGE_I??
	            avglc.update_mask_col(self.flags['flag_smallnum'], [avglc_index])
	        # else check sigmacut bounds and flag
	        else:
	            is_bad = False
	            if fluxstatparams['Ngood'] < self.Ngood_min:
	                is_bad = True
	            if fluxstatparams['Nclip'] > self.Nclip_max:
	                is_bad = True
	            if not(fluxstatparams['X2norm'] is None) and fluxstatparams['X2norm'] > self.x2_max:
	                is_bad = True
	            if is_bad:
	                lc.update_mask_col(self.flags['flag_badday'], range_i)
	                avglc.update_mask_col(self.flags['flag_badday'], [avglc_index])

	        mjd += self.mjd_bin_size
	    
	    # convert flux to magnitude and dflux to dmagnitude
	    avglc.pdastro.flux2mag('uJy','duJy','m','dm', zpt=23.9, upperlim_Nsigma=self.flux2mag_sigma_limit)

	    avglc = self.drop_extra_columns(avglc)

	    for col in ['Nclip','Ngood','Nexcluded','Mask']: 
	        avglc.pdastro.t[col] = avglc.pdastro.t[col].astype(np.int32)

		return lc, avglc

	def average_loop(self):
		args = self.define_args().parse_args()
		self.load_settings(args)

		for obj_index in range(0,len(args.tnsnames)):
			print(f'\nCOMMENCING AVERAGING LOOP FOR SN {args.tnsnames[obj_index]}')

			for filt in ['o','c']:
				print(f'\nFILTER SET: {filt}')
				lc = atlas_lc(tnsname=args.tnsnames[obj_index])
				lc.load(filt, self.input_dir)

				avglc = atlas_lc(tnsname=lc.tnsname,is_averaged=True,mjd_bin_size=self.mjd_bin_size)

				lc, avglc = average_lc(lc, avglc)

				lc.save(self.output_dir, filt=filt, overwrite=self.overwrite)
				avglc.save(self.output_dir, filt=filt, overwrite=self.overwrite)

if __name__ == "__main__":
	average_atlas_lc = average_atlas_lc()
	average_atlas_lc.average_loop()
