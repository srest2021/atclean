#!/usr/bin/env python
"""
@author: Sofia Rest
"""

import sys
from copy import deepcopys
import pandas as pd
import numpy as np

import sigmacut
from pdastro import pdastrostatsclass, AandB, AnotB
from light_curve import light_curve
from plot_lc import plot_lc

class cut_lc():
    def __init__(self):
    	# credentials
    	self.tns_api_key = None

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

		# chi-square cut
		self.chisquares = False
		self.chisquare_cut = None 
		self.stn_bound = None
		self.min_cut = None
		self.max_cut = None
		self.cut_step = None
		self.contam_lim = None
		self.loss_lim = None
		self.lim_to_prioritize = None

		# uncertainty cut
		self.uncertainties = False 
		self.uncertainty_cut = None

		# control light curve cut
		self.controls = False
		self.num_controls = None
		self.x2_max = None
		self.stn_max = None
		self.Nclip_max = None
		self.Ngood_min = None
    
    # define command line arguments
	def define_args(self, parser=None, usage=None, conflict_handler='resolve'):
		if parser is None:
			parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)
		
		parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')
		
		parser.add_argument('-x', '--chisquares', default=False, action='store_true', help='apply chi-square cut')
		parser.add_argument('-u', '--uncertainties', default=False, action='store_true', help='apply uncertainty cut')
		parser.add_argument('-c', '--controls', default=False, action='store_true', help='apply control light curve cut')
		
		parser.add_argument('-p', '--plot', default=False, action='store_true', help='plot each cut and save into PDF file')
		
		parser.add_argument('-f','--cfg_filename', default='atlaslc.ini', type=str, help='file name of ini file with settings for this class')
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

		self.chisquares = args.chisquares
		if self.chisquares:
			print(f'Chi-square cut: {self.chisquares}')
			if isdigit(cfg['Chi-square cut settings']['override_cut']):
				self.chisquare_cut = cfg['Chi-square cut settings']['override_cut']
				print(f'# Overriding dynamic chi-square cut with manual cut of x2 = {self.chisquare_cut}')
			else:
				self.stn_bound = cfg['Chi-square cut settings']['stn_bound']
				self.min_cut = cfg['Chi-square cut settings']['min_cut']
				self.max_cut = cfg['Chi-square cut settings']['max_cut']
				self.cut_step = cfg['Chi-square cut settings']['cut_step']
				self.contam_lim = cfg['Chi-square cut settings']['contamination_limit']
				self.loss_lim = cfg['Chi-square cut settings']['loss_limit']
				self.lim_to_prioritize = cfg['Chi-square cut settings']['limit_to_prioritize']
				if not(self.lim_to_prioritize is 'loss') and not(self.lim_to_prioritize is 'contamination'):
					raise RuntimeError('ERROR: Limit to prioritize (limit_to_prioritize in config file) must be set to \'contamination\' or \'loss\'!')
				print(f'# abs(flux/dflux) bound that determines a "good" measurement vs. "bad" measurement: {self.stn_bound}')
				print(f'# Cut range: [{self.min_cut}, {self.max_cut}], both ends inclusive, with step size {self.cut_step}')
				print(f'# Contamination percent limit: {self.contam_lim}')
				print(f'# Loss percent limit: {self.loss_lim}')
				print(f'# Limit to prioritize: {self.lim_to_prioritize}')

		self.uncertainties = args.uncertainties
		if self.uncertainties:
			print(f'Uncertainty cut: {self.uncertainties}')
			self.uncertainty_cut = cfg['Uncertainty cut settings']['cut']
			print(f'# Set to cut at dflux = {self.uncertainty_cut}')

		self.controls = args.controls
		if self.controls:
			print(f'Control light curve cut: {self.controls}')
			self.num_controls = cfg['Control light curve settings']['num_controls']
			self.x2_max = cfg['Control light curve settings']['x2_max']
			self.stn_max = cfg['Control light curve settings']['stn_max']
			self.Nclip_max = cfg['Control light curve settings']['Nclip_max']
			self.Ngood_min = cfg['Control light curve settings']['Ngood_min']
			print(f'# Bound for an epoch\'s maximum chi-square: {self.x2_max}')
			print(f'# Bound for an epoch\'s maximum abs(flux/dflux) ratio: {self.stn_max}')
			print(f'# Bound for an epoch\'s maximum number of clipped control measurements: {self.Nclip_max}')
			print(f'# Bound for an epoch\'s minimum number of good control measurements: {self.Ngood_min}')

		self.plot = args.plot

	# helper function for get_baseline_regions()
	def get_Ndays(self, SN_region_index):
		return 200 if SN_region_index == 2 else 40

	# get regions of a lc where no SN flux is present
	def get_baseline_regions(self, lc, Ndays_min):
		print('# Getting region indices around SN... ')

		baseline_ix = lc.get_baseline_ix()

		regions = {}
		regions['t0'] = lc.pdastro.ix_inrange(colnames=['MJD'], uplim=tchange1)
		regions['t1'] = lc.pdastro.ix_inrange(colnames=['MJD'], lowlim=tchange1, uplim=tchange2)
		regions['t2'] = lc.pdastro.ix_inrange(colnames=['MJD'], lowlim=tchange2)
		regions['b_t0'] = AandB(regions['t0'], baseline_ix)
		regions['b_t1'] = AandB(regions['t1'], baseline_ix)
		regions['b_t2'] = AandB(regions['t2'], baseline_ix)

		# find region SN starts in 
		SN_region_index = None
		for region_index in range(0,3):
			region_i = regions['t%d' % region_index]
			if lc.discdate >= lc.pdastro.t.loc[region_i[0],'MJD'] and lc.discdate <= lc.pdastro.t.loc[region_i[-1],'MJD']:
				SN_region_index = region_index
		if SN_region_index is None:
			raise RuntimeError('## ERROR: could not find region with SN discovery date!')
		else:
			print('## SN discovery date located in template region t%d' % SN_region_index)

		# for region with tail end of the SN, get last Ndays days and classify as baseline
		adjust_region_index = SN_region_index
		if adjust_region_index < 2 and len(regions['b_t%d'%adjust_region_index]) >= Ndays_min:
			adjust_region_index += 1
		if len(regions['b_t%d'%adjust_region_index]) < Ndays_min:
			print('## Getting baseline flux for template region t%d by obtaining last %d days of region... ' % (adjust_region_index, self.get_Ndays(adjust_region_index)))
			regions['b_t%d'%adjust_region_index] = lc.pdastro.ix_inrange(colnames=['MJD'],
																		lowlim=lc.pdastro.t.loc[regions['t%d'%adjust_region_index][-1],'MJD'] - self.get_Ndays(adjust_region_index),
																		uplim=lc.pdastro.t.loc[regions['t%d'%adjust_region_index][-1],'MJD'])
		if adjust_region_index < 1: regions['b_t1'] = regions['t1']
		if adjust_region_index < 2: regions['b_t2'] = regions['t2']

		for region_index in range(0,3):
			if len(regions['t%d'%region_index]) > 0:
				print('## TEMPLATE REGION t%d MJD RANGE: %0.2f - %0.2f' % (region_index, lc.pdastro.t.loc[regions['t%d'%region_index][0],'MJD'], lc.pdastro.t.loc[regions['t%d'%region_index][-1],'MJD']))
			else:
				print('## TEMPLATE REGION t%d MJD RANGE: not found' % region_index)
			if len(regions['b_t%d'%region_index]) > 0:
				print('## TEMPLATE REGION b_t%d BASELINE MJD RANGE: %0.2f - %0.2f' % (region_index, lc.pdastro.t.loc[regions['b_t%d'%region_index][0],'MJD'], lc.pdastro.t.loc[regions['b_t%d'%region_index][-1],'MJD']))
			else:
				print('## TEMPLATE REGION b_t%d BASELINE MJD RANGE: not found' % region_index)

		# check to make sure baseline flux is still consistent by getting median of first and last halves of affected region
		first_i = regions['b_t%d'%adjust_region_index][0]
		mid_i   = regions['b_t%d'%adjust_region_index][int(len(regions['b_t%d'%adjust_region_index])/2)]
		last_i  = regions['b_t%d'%adjust_region_index][-1]
		median1 = np.median(lc.pdastro.t.loc[lc.pdastro.ix_inrange(colnames=['MJD'], lowlim=lc.pdastro.t.loc[first_i,'MJD'], uplim=lc.pdastro.t.loc[mid_i,'MJD']), 'uJy'])
		median2 = np.median(lc.pdastro.t.loc[lc.pdastro.ix_inrange(colnames=['MJD'], lowlim=lc.pdastro.t.loc[mid_i+1,'MJD'], uplim=lc.pdastro.t.loc[last_i,'MJD']), 'uJy'])
		print(f'## Checking that baseline flux is consistent throughout adjusted region...\n## Median of first half: {median1:%0.2f}\n## Median of second half: {median2:%0.2f}')

		lc.corrected_baseline_ix = np.concatenate([regions['b_t0'], regions['b_t1'], regions['b_t2']])
		lc.during_sn_ix = AnotB(lc.pdastro.getindices(), lc.corrected_baseline_ix)
		
		return regions, lc

	# correct for atlas template changes at mjd=58417,58882 
	# more info here: https://fallingstar-data.com/forcedphot/faq/
	def correct_for_template(self, lc):
		print('Correcting for potential flux in template due to template changes at MJD=58417,58882...')
		
		# automatically define baseline regions according to discovery date
		regions, lc = self.get_baseline_regions(lc, Ndays_min=6)

		# get indices of measurements with x2<=5 so that when getting median, use these indices if possible
		b_goodx2_i = lc.pdastro.ix_inrange(colnames=['chi/N'], uplim=5, indices=lc.corrected_baseline_ix)

		# for each region, adjust for template change by subtracting median of that region's baseline flux
		for region_index in range(0,3):
			region_i = regions['b_t%d'%region_index]
			if len(region_i) > 0:
				print(f'# Adjusting for template change in region b_t{region_index:%d} from {lc.pdastro.t.loc[region_i[0],'MJD']:%0.2f}-{lc.pdastro.t.loc[region_i[-1],'MJD']:%0.2f}...')
				print(f'## Baseline median before: {np.median(lc.pdastro.t.loc[region_i,'uJy'])}')
				if len(AandB(region_i,b_goodx2_i)) > 0:
					median = np.median(lc.pdastro.t.loc[AandB(region_i,b_goodx2_i),'uJy'])
				else:
					median = np.median(lc.pdastro.t.loc[region_i,'uJy'])
				print(f'## Subtracting median {median:0.1f} uJy of baseline flux with chi-square ≤ 5 from light curve flux due to potential flux in the template...')
				lc.pdastro.t.loc[regions['t%d'%region_index],'uJy'] -= median
				print(f'## Baseline median now: {np.median(lc.pdastro.t.loc[region_i,'uJy'])}')
			else:
				print(f'# No baseline region for region b_t{region_index}, skipping...')

		return lc

	# for a range of chi-square cuts, determine contamination, loss, and other pecentages
	def get_limcuts_table(self, lc, indices=None):
		limcuts = pdastrostatsclass(columns=['PSF Chi-Square Cut', 'N', 'Ngood', 'Nbad', 'Nkept', 'Ncut', 'Ngood,kept', 'Ngood,cut', 'Nbad,kept', 'Nbad,cut',
                                          	  'Pgood,kept', 'Pgood,cut', 'Pbad,kept', 'Pbad,cut', 'Ngood,kept/Ngood', 'Ploss', 'Pcontamination',
                                              'Nbad,cut 3<stn<=5', 'Nbad,cut 5<stn<=10', 'Nbad,cut 10<stn', 'Nbad,kept 3<stn<=5', 'Nbad,kept 5<stn<=10', 'Nbad,kept 10<stn'])

		if indices is None:
			indices = lc.corrected_baseline_ix
		
		# define good and bad baseline measurement indices according to abs(uJy/duJy) bound
	    b_good_i = lc.pdastro.ix_inrange(colnames=['uJy/duJy'], lowlim=-self.stn_bound, uplim=self.stn_bound, indices=indices)
	    b_bad_i = AnotB(indices, b_good_i)

	    # for different x2 cuts decreasing from 50
	    for cut in range(self.min_cut, self.max_cut+1, self.cut_step):
	        # define kept baseline measurement indices according to current chi-square cut
	        b_kept_i = lc.pdastro.ix_inrange(colnames=['chi/N'], uplim=cut, indices=indices)
	        b_cut_i = AnotB(indices, b_kept_i)

	        # construct new row of limcuts table
	        df = pd.DataFrame([[cut, len(indices), # N
                            len(b_good_i), # Ngood
                            len(b_bad_i), # Nbad
                            len(b_kept_i), # Nkept
                            len(b_cut_i), # Ncut
                            len(AandB(b_good_i,b_kept_i)), # Ngood,kept
                            len(AandB(b_good_i,b_cut_i)), # Ngood,cut
                            len(AandB(b_bad_i,b_kept_i)), # Nbad,kept
                            len(AandB(b_bad_i,b_cut_i)), # Nbad,cut
                            100*len(AandB(b_good_i,b_kept_i))/len(indices), # Ngood,kept/Nbaseline
                            100*len(AandB(b_good_i,b_cut_i))/len(indices), # Ngood,cut/Nbaseline 
                            100*len(AandB(b_bad_i,b_kept_i))/len(indices), # Nbad,kept/Nbaseline
                            100*len(AandB(b_bad_i,b_cut_i))/len(indices), # Nbad,cut/Nbaseline
                            100*len(AandB(b_good_i,b_kept_i))/len(b_good_i), # Ngood,kept/Ngood
                            100*len(AandB(b_good_i,b_cut_i))/len(b_good_i), # Ngood,cut/Ngood = Loss
                            100*len(AandB(b_bad_i,b_kept_i))/len(b_kept_i), # Nbad,kept/Nkept = Contamination
                            len(AandB(AandB(b_bad_i,b_cut_i), lc.pdastro.ix_inrange(colnames=['uJy/duJy'], lowlim=-3, uplim=5, exclude_lowlim=True))), # Nbad,cut 3<stn<=5
                            len(AandB(AandB(b_bad_i,b_cut_i), lc.pdastro.ix_inrange(colnames=['uJy/duJy'], lowlim=5, uplim=10, exclude_lowlim=True))), # Nbad,cut 5<stn<=10
                            len(AandB(AandB(b_bad_i,b_cut_i), lc.pdastro.ix_inrange(colnames=['uJy/duJy'], lowlim=10, exclude_lowlim=True))), # Nbad,cut 10<stn 
                            len(AandB(AandB(b_bad_i,b_kept_i), lc.pdastro.ix_inrange(colnames=['uJy/duJy'], lowlim=-3, uplim=5, exclude_lowlim=True))), # Nbad,kept 3<stn<=5
                            len(AandB(AandB(b_bad_i,b_kept_i), lc.pdastro.ix_inrange(colnames=['uJy/duJy'], lowlim=5, uplim=10, exclude_lowlim=True))), # Nbad,kept 5<stn<=10
                            len(AandB(AandB(b_bad_i,b_kept_i), lc.pdastro.ix_inrange(colnames=['uJy/duJy'], lowlim=10, exclude_lowlim=True))), # Nbad,kept 10<stn 
                            ]], columns=['PSF Chi-Square Cut', 'N', 'Ngood', 'Nbad', 'Nkept', 'Ncut', 'Ngood,kept', 'Ngood,cut', 'Nbad,kept', 'Nbad,cut',
                                         'Pgood,kept', 'Pgood,cut', 'Pbad,kept', 'Pbad,cut', 'Ngood,kept/Ngood', 'Ploss', 'Pcontamination',
                                         'Nbad,cut 3<stn<=5', 'Nbad,cut 5<stn<=10', 'Nbad,cut 10<stn', 'Nbad,kept 3<stn<=5', 'Nbad,kept 5<stn<=10', 'Nbad,kept 10<stn'])
        	limcuts.t = pd.concat([limcuts.t,df],ignore_index=True)
		return limcuts 

	# get contamination and loss plus other percentages and information for a certain chi-square cut 
	def get_limcuts_data(self, lc, name, cut, case):
		indices = lc.corrected_baseline_ix
	  
	    b_good_i = lc.pdastro.ix_inrange(colnames=['uJy/duJy'],lowlim=-stn_bound,uplim=stn_bound,indices=indices)
	    b_bad_i = AnotB(indices, b_good_i)
	    b_kept_i = lc.pdastro.ix_inrange(colnames=['chi/N'],uplim=cut,indices=indices)
	    b_cut_i = AnotB(indices, b_kept_i)

	    data = {}
	    data['name'] = name
	    data['cut'] = cut
	    data['case'] = case
	    data['Ngood'] = len(b_good_i)
	    data['Nbad'] = len(b_bad_i)
	    data['Nkept'] = len(b_kept_i)
	    data['Ncut'] = len(b_cut_i)
	    data['Ngood,kept'] = len(AandB(b_good_i,b_kept_i))
	    data['Ngood,cut'] = len(AandB(b_good_i,b_cut_i))
	    data['Nbad,kept'] = len(AandB(b_bad_i,b_kept_i))
	    data['Nbad,cut'] = len(AandB(b_bad_i,b_cut_i))
	    data['Pgood,kept'] = 100*len(AandB(b_good_i,b_kept_i))/len(indices)
	    data['Pgood,cut'] = 100*len(AandB(b_good_i,b_cut_i))/len(indices)
	    data['Pbad,kept'] = 100*len(AandB(b_bad_i,b_kept_i))/len(indices)
	    data['Pbad,cut'] = 100*len(AandB(b_bad_i,b_cut_i))/len(indices)
	    data['Ngood,kept/Ngood'] = 100*len(AandB(b_good_i,b_kept_i))/len(b_good_i)
	    data['Ploss'] = 100*len(AandB(b_good_i,b_cut_i))/len(b_good_i)
	    data['Pcontamination'] = 100*len(AandB(b_bad_i,b_kept_i))/len(b_kept_i)

	# get the optimal contamination and loss cuts according to percentages in limcuts table
	def get_limcuts(self, limcuts):
		contam_cut = None
	    loss_cut = None
	    contam_case = None
	    loss_case = None

	    sortby_loss = limcuts.t.iloc[(limcuts.t['Ploss']).argsort()].reset_index()
	    min_loss = sortby_loss.loc[0,'Ploss']
	    max_loss = sortby_loss.loc[len(sortby_loss)-1,'Ploss']
	    # if all loss below lim, loss_cut is min cut
	    if min_loss < loss_lim and max_loss < loss_lim:
	        loss_case = 'below lim'
	        loss_cut = limcuts.t.loc[0,'PSF Chi-Square Cut']
	    else:
	        # else if all loss above lim, loss_cut is min cut with min% loss
	        if min_loss > loss_lim and max_loss > loss_lim:
	            loss_case = 'above lim'
	            a = np.where(limcuts.t['Ploss'] == min_loss)[0]
	            b = limcuts.t.iloc[a]
	            c = b.iloc[(b['PSF Chi-Square Cut']).argsort()].reset_index()
	            loss_cut = c.loc[0,'PSF Chi-Square Cut']
	        # else if loss crosses lim at some point, loss_cut is min cut with max% loss <= loss_lim
	        else:
	            loss_case = 'crosses lim'
	            valid_cuts = sortby_loss[sortby_loss['Ploss'] <= loss_lim]
	            a = np.where(limcuts.t['Ploss'] == valid_cuts.loc[len(valid_cuts)-1,'Ploss'])[0]
	            # sort by cuts
	            b = limcuts.t.iloc[a]
	            c = b.iloc[(b['PSF Chi-Square Cut']).argsort()].reset_index()
	            # get midpoint of loss1 and loss2 (two points on either side of lim)
	            loss1_i = np.where(limcuts.t['PSF Chi-Square Cut'] == c.loc[0,'PSF Chi-Square Cut'])[0][0]
	            if limcuts.t.loc[loss1_i,'Ploss'] == loss_lim:
	                loss_cut = limcuts.t.loc[loss1_i,'PSF Chi-Square Cut']
	            else:
	                loss2_i = loss1_i - 1
	                x = np.array([limcuts.t.loc[loss1_i,'PSF Chi-Square Cut'], limcuts.t.loc[loss2_i,'PSF Chi-Square Cut']])
	                contam_y = np.array([limcuts.t.loc[loss1_i,'Pcontamination'], limcuts.t.loc[loss2_i,'Pcontamination']])
	                loss_y = np.array([limcuts.t.loc[loss1_i,'Ploss'], limcuts.t.loc[loss2_i,'Ploss']])
	                contam_line = np.polyfit(x,contam_y,1)
	                loss_line = np.polyfit(x,loss_y,1)
	                loss_cut = (loss_lim-loss_line[1])/loss_line[0]

	    sortby_contam = limcuts.t.iloc[(limcuts.t['Pcontamination']).argsort()].reset_index()
	    min_contam = sortby_contam.loc[0,'Pcontamination']
	    max_contam = sortby_contam.loc[len(sortby_contam)-1,'Pcontamination']
	    # if all contam below lim, contam_cut is max cut
	    if min_contam < contam_lim and max_contam < contam_lim:
	        contam_case = 'below lim'
	        contam_cut = limcuts.t.loc[len(limcuts.t)-1,'PSF Chi-Square Cut']
	    else:
	        # else if all contam above lim, contam_cut is max cut with min% contam
	        if min_contam > contam_lim and max_contam > contam_lim:
	            contam_case = 'above lim'
	            a = np.where(limcuts.t['Pcontamination'] == min_contam)[0]
	            b = limcuts.t.iloc[a]
	            c = b.iloc[(b['PSF Chi-Square Cut']).argsort()].reset_index()
	            contam_cut = c.loc[len(c)-1,'PSF Chi-Square Cut']
	        # else if contam crosses lim at some point, contam_cut is max cut with max% contam <= contam_lim
	        else:
	            contam_case = 'crosses lim'
	            valid_cuts = sortby_contam[sortby_contam['Pcontamination'] <= contam_lim]
	            a = np.where(limcuts.t['Pcontamination'] == valid_cuts.loc[len(valid_cuts)-1,'Pcontamination'])[0]
	            # sort by cuts
	            b = limcuts.t.iloc[a]
	            c = b.iloc[(b['PSF Chi-Square Cut']).argsort()].reset_index()
	            # get midpoint of contam1 and contam2 (two points on either side of lim)
	            contam1_i = np.where(limcuts.t['PSF Chi-Square Cut'] == c.loc[len(c)-1,'PSF Chi-Square Cut'])[0][0]
	            if limcuts.t.loc[contam1_i,'Pcontamination'] == contam_lim:
	                contam_cut = limcuts.t.loc[contam1_i,'PSF Chi-Square Cut']
	            else:
	                contam2_i = contam1_i + 1
	                x = np.array([limcuts.t.loc[contam1_i,'PSF Chi-Square Cut'], limcuts.t.loc[contam2_i,'PSF Chi-Square Cut']])
	                contam_y = np.array([limcuts.t.loc[contam1_i,'Pcontamination'], limcuts.t.loc[contam2_i,'Pcontamination']])
	                loss_y = np.array([limcuts.t.loc[contam1_i,'Ploss'], limcuts.t.loc[contam2_i,'Ploss']])
	                contam_line = np.polyfit(x,contam_y,1)
	                loss_line = np.polyfit(x,loss_y,1)
	                contam_cut = (contam_lim-contam_line[1])/contam_line[0]

	    loss_cut_data = self.get_limcuts_data(lc, 'loss_cut', loss_cut, loss_case)
	    contam_cut_data = self.get_limcuts_data(lc, 'contam_cut', contam_cut, contam_case)

	    return loss_cut_data, contam_cut_data

	# choose between contamination cut and loss cut
	def get_final_chisquare_cut(self, contam_cut, loss_cut, contam_case, loss_case):
		# case 1 and 1: final_cut = 3
	    # case 1 and 2: take limit of case 2
	    # case 1 and 3: take limit of case 3
	    # case 2 and 2: print lims don't work
	    # case 2 and 3: get_final_chisquare_cut
	    # case 3 and 3: get_final_chisquare_cut

	    case1 = loss_case == 'below lim' or contam_case == 'below lim'
	    case2 = loss_case == 'above lim' or contam_case == 'above lim'
	    case3 = loss_case == 'crosses lim' or contam_case == 'crosses lim'

	    final_cut = None
	    if case1 and not case2 and not case3: # 1 and 1
	        print('Valid chi-square cut range from %0.2f to %0.2f! Setting to 3...' % (loss_lim_cut, contam_lim_cut))
	        final_cut = cut_start
	    elif case1: # 1
	        if case2: # and 2
	            if loss_case == 'above lim':
	                print('WARNING: contam_lim_cut <= %0.2f falls below limit %0.2f%%, but loss_lim_cut >= %0.2f falls above limit %0.2f%%! Setting to %0.2f...' % (contam_lim_cut, contam_lim, loss_lim_cut, loss_lim, loss_lim_cut))
	                final_cut = loss_lim_cut
	            else:
	                print('WARNING: loss_lim_cut <= %0.2f falls below limit %0.2f%%, but contam_lim_cut >= %0.2f falls above limit %0.2f%%! Setting to %0.2f...' % (loss_lim_cut, loss_lim, contam_lim_cut, contam_lim, contam_lim_cut))
	                final_cut = contam_lim_cut
	        else: # and 3
	            if loss_case == 'crosses lim':
	                print('Contam_lim_cut <= %0.2f falls below limit %0.2f%% and loss_lim_cut >= %0.2f crosses limit %0.2f%%, setting to %0.2f...' % (contam_lim_cut, contam_lim, loss_lim_cut, loss_lim, loss_lim_cut))
	                final_cut = loss_lim_cut
	            else:
	                print('Loss_lim_cut <= %0.2f falls below limit %0.2f%% and contam_lim_cut >= %0.2f crosses limit %0.2f%%, setting to %0.2f...' % (loss_lim_cut, loss_lim, contam_lim_cut, contam_lim, contam_lim_cut))
	                final_cut = contam_lim_cut
	    elif case2 and not case3: # 2 and 2
	        print('ERROR: chi-square loss_lim_cut >= %0.2f and contam_lim_cut <= %0.2f both fall above limits %0.2f%% and %0.2f%%! Try setting less strict limits. Setting final cut to nan.' % (loss_lim_cut, contam_lim_cut, loss_lim, contam_lim))
	        final_cut = np.nan
	    else: # 2 and 3 or 3 and 3
	        if loss_lim_cut > contam_lim_cut:
	            print('WARNING: chi-square loss_lim_cut >= %0.2f and contam_lim_cut <= %0.2f do not overlap! ' % (loss_lim_cut, contam_lim_cut))
	            if lim_to_prioritize == 'contam_lim':
	                print('Prioritizing %s and setting to %0.2f...' % (lim_to_prioritize, contam_lim_cut))
	                final_cut = contam_lim_cut
	            else:
	                print('Prioritizing %s and setting to %0.2f... ' % (lim_to_prioritize, loss_lim_cut))
	                final_cut = loss_lim_cut
	        else:
	            print('Valid chi-square cut range from %0.2f to %0.2f! ' % (loss_lim_cut, contam_lim_cut))
	            if lim_to_prioritize == 'contam_lim':
	                print('Prioritizing %s and setting to %0.2f... ' % (lim_to_prioritize, loss_lim_cut))
	                final_cut = loss_lim_cut
	            else:
	                print('Prioritizing %s and setting to %0.2f... ' % (lim_to_prioritize, contam_lim_cut))
	                final_cut = contam_lim_cut
	    return final_cut

	# apply chi-square cut to SN light curve and update mask column with flag
	def apply_chisquare_cut(self, lc):
		print('Now applying dynamic chi-square cut...')
		
		limcuts = self.get_limcuts_table(lc)

		loss_cut_data, contam_cut_data = self.get_limcuts(limcuts)

	    print(f'# Contamination cut according to given contam_limit, with {contam_cut_data['Pcontamination']:0.2f}% contamination and {contam_cut_data['Ploss']:0.2f}% loss: {contam_cut_data['cut']:0.2f}')
		if contam_cut_data['case'] == 'above lim':
		    print(f'## WARNING: Contamination cut not possible with contamination <= contam_lim {self.contam_lim:0.2f}%!')
		print(f'# Loss cut according to given loss_limit, with {loss_cut_data['Pcontamination']:0.2f}% contamination and {loss_cut_data['Ploss']:0.2f}% loss: {loss_cut_data['cut']:0.2f}')
		if loss_cut_data['case'] == 'above lim':
		    print(f'## WARNING: Loss cut not possible with loss <= loss_lim {self.loss_lim:0.2f}%!')

		self.chisquare_cut = self.get_final_chisquare_cut(contam_cut_data['cut'], loss_cut_data['cut'], contam_cut_data['case'], loss_cut_data['case'])

		if np.isnan(self.chisquare_cut):
		    raise RuntimeError('\n# ERROR: Final suggested chi-square cut could not be determined according to given contamination and loss limits. We suggest resetting your limits in atlaslc.ini.')
		else:
		    if self.chisquare_cut == contam_cut:
		    	Pcontamination = contam_cut_data['Pcontamination']
		    	Ploss = contam_cut_data['Ploss']
		    else:
		    	Pcontamination = loss_cut_data['Pcontamination']
		    	Ploss = loss_cut_data['Ploss']
		    print(f'# Final suggested chi-square cut is {self.chisquare_cut:0.2f}, with {Pcontamination:0.2f}% contamination and {Ploss:0.2f}% loss.')
		   	if (Pcontamination > contam_lim):
		        print(f'## WARNING: Final cut\'s contamination {Pcontamination:0.2f}% exceeds {contam_lim:0.2f}%!')
		    if (Ploss > loss_lim):
		        print(f'## WARNING: Final cut\'s loss {Ploss:0.2f}% exceeds loss_lim {loss_lim:0.2f}%!')

		# remove old mask column
		if 'Mask' in lc.pdastro.t.columns: 
		    lc.pdastro.t.drop(columns=['Mask'],inplace=True)

		# create new mask column and update it with final chi-square cut
		lc.pdastro.t['Mask'] = 0
		cut_ix = lc.pdastro.ix_inrange(colnames=['chi/N'], lowlim=self.chisquare_cut, exclude_lowlim=True)
		lc.update_mask_col(flags['flag_chisquare'], cut_ix)
		print(f'# Total % of data cut: {len(cut_ix)/len(lc.pdastro.getindices())}%')

		return lc

	# apply chi-square cut to SN light curve and update mask column with flag
	def apply_uncertainty_cut(self, lc):
		print('Now applying uncertainty cut...')

		cut_ix = lc.pdastro.ix_inrange(colnames=['duJy'], lowlim=self.uncertainty_cut, exclude_lowlim=True)
		lc.update_mask_col(flags['flag_uncertainty'], cut_ix)
		print(f'# Total % of data cut: {len(cut_ix)/len(lc.pdastro.getindices())}%')

		return lc

	def cut_loop(self):
		args = self.define_args().parse_args()
		self.load_settings(args)

		for obj_index in range(0,len(args.tnsnames)):
			print(f'\nCommencing cut loop for SN {lc.tnsname}')

			if args.plot:
				plot_lc = plot_lc(output_dir=self.output_dir, flags=self.flags)

			for filt in ['o','c']:
				print(f'Filter set: {filt}')
				lc = light_curve(tnsname=args.tnsnames[obj_index])
				lc.load(filt, self.input_dir, num_controls=self.num_controls)
				lc.get_tns_data(self.tns_api_key)
				lc = self.correct_for_template(lc)

				if self.chisquares and self.chisquare_cut is None:
					lc = self.apply_chisquare_cut(lc)

				if self.uncertainties:
					lc = self.apply_uncertainty_cut()


				


				lc.save(self.output_dir, filt=filt, overwrite=self.overwrite)

				if args.plot:
					plot_lc.set(lc=lc, filt=filt)
					plot_lc.plot_lc()
					plot_lc.plot_chisquare_cut()
					plot_lc.plot_uncertainty_cut()
					plot_lc.plot_controls_cut()

if __name__ == "__main__":
	cut_lc = cut_lc()
	cut_lc.cut_loop()