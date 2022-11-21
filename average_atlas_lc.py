#!/usr/bin/env python
"""
@author: Sofia Rest
"""

import sys, argparse, configparser, re, os
from copy import deepcopy
import pandas as pd
import numpy as np
import PyPDF2 

import sigmacut
from pdastro import pdastrostatsclass, AandB, AnotB
from atlas_lc import atlas_lc
from plot_atlas_lc import plot_atlas_lc

class average_atlas_lc():
	def __init__(self):
		# credentials
		self.tns_api_key = None

		# input/output
		self.input_dir = None
		self.output_dir = None
		self.overwrite = True

		# flags
		self.flags = {'chisquare':0x1,

					  'uncertainty':0x2,

					  'controls_bad':0x400000,
					  'controls_questionable':0x80000,

					  'avg_badday':0x800000,
					  'avg_ixclip':0x1000,
					  'avg_smallnum':0x2000}

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
		parser.add_argument('-a','--tns_api_key', type=str, help='api key to access TNS')

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

		self.tns_api_key = cfg['TNS credentials']['api_key'] if args.tns_api_key is None else args.tns_api_key

		self.input_dir = cfg['Input/output settings']['input_dir']
		print(f'Light curve .txt files input directory: {self.input_dir}')
		self.output_dir = cfg['Input/output settings']['output_dir']
		print(f'Light curve .txt files output directory: {self.output_dir}')
		self.overwrite = not args.dont_overwrite
		print(f'Overwrite existing light curve files: {self.overwrite}')
		self.flux2mag_sigmalimit = int(cfg['Input/output settings']['flux2mag_sigmalimit'])
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

	# helper function for set_corrected_baseline_ix()
	def get_Ndays(self, SN_region_index):
		return 200 if SN_region_index == 2 else 40

	def set_corrected_baseline_ix(self, avglc):
		print('\nGetting region indices around SN... ')
		Ndays_min = 6

		baseline_ix = avglc.get_baseline_ix()
		tchange1 = 58417
		tchange2 = 58882

		regions = {}
		regions['t0'] = avglc.pdastro.ix_inrange(colnames=['MJD'], uplim=tchange1)
		regions['t1'] = avglc.pdastro.ix_inrange(colnames=['MJD'], lowlim=tchange1, uplim=tchange2)
		regions['t2'] = avglc.pdastro.ix_inrange(colnames=['MJD'], lowlim=tchange2)
		regions['b_t0'] = AandB(regions['t0'], baseline_ix)
		regions['b_t1'] = AandB(regions['t1'], baseline_ix)
		regions['b_t2'] = AandB(regions['t2'], baseline_ix)

		# find region SN starts in 
		SN_region_index = None
		if avglc.discdate <= tchange1:
			SN_region_index = 0
		elif avglc.discdate > tchange1 and avglc.discdate <= tchange2:
			SN_region_index = 1
		elif avglc.discdate > tchange2:
			SN_region_index = 2
		if SN_region_index is None:
			raise RuntimeError('# ERROR: Something went wrong--could not find region with SN discovery date!')
		else:
			print('# SN discovery date located in template region t%d' % SN_region_index)

		# for region with tail end of the SN, get last Ndays days and classify as baseline
		adjust_region_index = SN_region_index
		if adjust_region_index < 2 and len(regions['b_t%d'%adjust_region_index]) >= Ndays_min:
			adjust_region_index += 1
		if len(regions['b_t%d'%adjust_region_index]) < Ndays_min:
			print('# Getting baseline flux for template region t%d by obtaining last %d days of region... ' % (adjust_region_index, self.get_Ndays(adjust_region_index)))
			regions['b_t%d'%adjust_region_index] = avglc.pdastro.ix_inrange(colnames=['MJD'],
																			lowlim=avglc.pdastro.t.loc[regions['t%d'%adjust_region_index][-1],'MJD'] - self.get_Ndays(adjust_region_index),
																			uplim=avglc.pdastro.t.loc[regions['t%d'%adjust_region_index][-1],'MJD'])
		if adjust_region_index < 1: regions['b_t1'] = regions['t1']
		if adjust_region_index < 2: regions['b_t2'] = regions['t2']

		avglc.corrected_baseline_ix = np.concatenate([regions['b_t0'], regions['b_t1'], regions['b_t2']])
		avglc.during_sn_ix = AnotB(avglc.pdastro.getindices(), avglc.corrected_baseline_ix)

		return avglc

	# add averaged light curve plots to original plot file
	# TO DO: does last plot in the original plot file need to be redone (all cuts plot)?
	def plot_averaged_lc(self, args, avglc, filt):
		plot = plot_atlas_lc(tnsname=avglc.tnsname, output_dir=self.output_dir, args=args, add2filename='avg', flags=self.flags)
		plot.set(lc=avglc, filt=filt)
		plot.plot_averaged_lc()
		plot.save()

		print('Merging original light curve plot PDF with averaged light curve plot PDF...')

		og_filename = f'{self.output_dir}/{avglc.tnsname}/{avglc.tnsname}_plots.pdf'
		avg_filename = f'{self.output_dir}/{avglc.tnsname}/{avglc.tnsname}_plots_avg.pdf'

		merger = PyPDF2.PdfMerger()
		try:
			merger.append(PyPDF2.PdfFileReader(og_filename, 'rb'))
		except Exception as e:
			print('No original plots...') # delete me
			pass
		merger.append(PyPDF2.PdfFileReader(avg_filename, 'rb'))
		merger.write(og_filename)
		merger.close()

		print('Removing averaged light curve plot PDF...')
		os.remove(avg_filename)

	def average_lc(self, lc, avglc):
		print(f'Averaging SN light curve...')

		mjd = int(np.amin(lc.pdastro.t['MJD']))
		mjd_max = int(np.amax(lc.pdastro.t['MJD']))+1

		good_i = lc.pdastro.ix_unmasked('Mask', maskval=self.flags['chisquare']|self.flags['uncertainty']|self.flags['controls_bad'])

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
					avglc.update_mask_col(self.flags['avg_badday'], [avglc_index])
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
				lc.update_mask_col(self.flags['avg_badday'], range_i)
				avglc.update_mask_col(self.flags['avg_badday'], [avglc_index])

				mjd += self.mjd_bin_size
				continue
			
			# average good measurements
			lc.pdastro.calcaverage_sigmacutloop('uJy', noisecol='duJy', indices=range_good_i, Nsigma=3.0, median_firstiteration=True)
			fluxstatparams = deepcopy(lc.pdastro.statparams)

			if fluxstatparams['mean'] is None or len(fluxstatparams['ix_good']) < 1:
				lc.update_mask_col(self.flags['avg_badday'], range_i)
				avglc.update_mask_col(self.flags['avg_badday'], [avglc_index])
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
				lc.update_mask_col(self.flags['avg_ixclip'], fluxstatparams['ix_clip'])
			
			# if small number within this bin, flag measurements
			if len(range_good_i) < 3:
				lc.update_mask_col(self.flags['avg_smallnum'], range_good_i) # TO DO: CHANGE TO RANGE_I??
				avglc.update_mask_col(self.flags['avg_smallnum'], [avglc_index])
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
					lc.update_mask_col(self.flags['avg_badday'], range_i)
					avglc.update_mask_col(self.flags['avg_badday'], [avglc_index])

			mjd += self.mjd_bin_size
		
		# convert flux to magnitude and dflux to dmagnitude
		avglc.pdastro.flux2mag('uJy','duJy','m','dm', zpt=23.9, upperlim_Nsigma=self.flux2mag_sigmalimit)

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
				lc.get_tns_data(self.tns_api_key)

				avglc = atlas_lc(tnsname=lc.tnsname,is_averaged=True,mjd_bin_size=self.mjd_bin_size)
				avglc.discdate = lc.discdate

				lc, avglc = self.average_lc(lc, avglc)

				avglc = self.set_corrected_baseline_ix(avglc)
				#print(avglc.corrected_baseline_ix)
				#print(avglc.during_sn_ix)

				lc.save(self.output_dir, filt=filt, overwrite=self.overwrite)
				avglc.save(self.output_dir, filt=filt, overwrite=self.overwrite)

				if args.plot:
					self.plot_averaged_lc(args, avglc, filt)

if __name__ == "__main__":
	average_atlas_lc = average_atlas_lc()
	average_atlas_lc.average_loop()
