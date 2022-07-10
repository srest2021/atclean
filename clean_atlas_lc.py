#!/usr/bin/env python
"""
@author: Sofia Rest
"""

import sys, argparse, configparser, re, os
from copy import deepcopy
import pandas as pd
import numpy as np

import sigmacut
from pdastro import pdastrostatsclass, AandB, AnotB
from atlas_lc import atlas_lc
from plot_atlas_lc import plot_atlas_lc

class clean_atlas_lc():
	def __init__(self):
		# credentials
		self.tns_api_key = None

		# input/output
		self.output_dir = None
		self.snlist_filename = None
		self.snlist = None
		self.overwrite = True

		# flags for each cut
		self.flags = {'chisquare':0x1,

					  'uncertainty':0x2,

					  'controls_bad':0x400000,
					  'controls_questionable':0x80000,
					  'controls_x2':0x100,
					  'controls_stn':0x200,
					  'controls_Nclip':0x400,
					  'controls_Ngood':0x800,

					  'avg_badday':0x800000,
					  'avg_ixclip':0x1000,
					  'avg_smallnum':0x2000}

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
		self.c_x2_max = None
		self.stn_max = None
		self.c_Nclip_max = None
		self.c_Ngood_min = None

		# averaging
		self.mjd_bin_size = None
		self.keep_empty_bins = False
		self.flux2mag_sigmalimit = None
		self.g_Nclip_max = None
		self.g_Ngood_min = None 
		self.g_x2_max = None
	
	# define command line arguments
	def define_args(self, parser=None, usage=None, conflict_handler='resolve'):
		if parser is None:
			parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)
		
		parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')
		
		parser.add_argument('-x', '--chisquares', default=False, action='store_true', help='apply chi-square cut')
		parser.add_argument('-u', '--uncertainties', default=False, action='store_true', help='apply uncertainty cut')
		parser.add_argument('-c', '--controls', default=False, action='store_true', help='apply control light curve cut')
		
		parser.add_argument('-g', '--average', default=False, action='store_true', help='average light curves and cut bad days')
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
		self.output_dir = cfg['Input/output settings']['output_dir']
		print(f'Light curve .txt files output directory: {self.output_dir}')

		# attempt to load snlist.txt; if does not exist, create new snlist table
		self.snlist_filename = cfg['Input/output settings']['snlist_filename']
		if os.path.exists(self.snlist_filename):
			self.snlist = pdastrostatsclass()
			self.snlist.load_spacesep(f'{self.output_dir}/{self.snlist_filename}', delim_whitespace=True)
		else:
			self.snlist = pdastrostatsclass(columns=['tnsname', 'ra', 'dec', 'discovery_date', 'closebright_ra', 'closebright_dec'])

		self.overwrite = not args.dont_overwrite
		print(f'Overwrite existing light curve files: {self.overwrite}')
		print(f'Plotting: {args.plot}')

		self.chisquares = args.chisquares
		if self.chisquares:
			print(f'\nChi-square cut: {self.chisquares}')
			if cfg['Chi-square cut settings']['override_cut'].isdigit():
				self.chisquare_cut = cfg['Chi-square cut settings']['override_cut']
				print(f'# Overriding dynamic chi-square cut with manual cut of x2 = {self.chisquare_cut}')
			else:
				self.stn_bound = float(cfg['Chi-square cut settings']['stn_bound'])
				self.min_cut = int(cfg['Chi-square cut settings']['min_cut'])
				self.max_cut = int(cfg['Chi-square cut settings']['max_cut'])
				self.cut_step = int(cfg['Chi-square cut settings']['cut_step'])
				self.contam_lim = float(cfg['Chi-square cut settings']['contamination_limit'])
				self.loss_lim = float(cfg['Chi-square cut settings']['loss_limit'])
				self.lim_to_prioritize = cfg['Chi-square cut settings']['limit_to_prioritize']
				if not(self.lim_to_prioritize == 'loss') and not(self.lim_to_prioritize == 'contamination'):
					raise RuntimeError(f'ERROR: Limit to prioritize (limit_to_prioritize in config file) must be set to \'contamination\' or \'loss\' but currently set to {self.lim_to_prioritize}!')
				print(f'# abs(flux/dflux) bound that determines a "good" measurement vs. "bad" measurement: {self.stn_bound}')
				print(f'# Cut range: [{self.min_cut}, {self.max_cut}], both ends inclusive, with step size {self.cut_step}')
				print(f'# Contamination percent limit: {self.contam_lim}')
				print(f'# Loss percent limit: {self.loss_lim}')
				print(f'# Limit to prioritize: {self.lim_to_prioritize}')

		self.uncertainties = args.uncertainties
		if self.uncertainties:
			print(f'\nUncertainty cut: {self.uncertainties}')
			self.uncertainty_cut = float(cfg['Uncertainty cut settings']['cut'])
			print(f'# Set to cut at dflux = {self.uncertainty_cut}')

		self.controls = args.controls
		if self.controls:
			print(f'\nControl light curve cut: {self.controls}')
			self.num_controls = int(cfg['Control light curve settings']['num_controls'])
			self.c_x2_max = float(cfg['Control light curve cut settings']['x2_max'])
			self.stn_max = float(cfg['Control light curve cut settings']['stn_max'])
			self.c_Nclip_max = int(cfg['Control light curve cut settings']['Nclip_max'])
			self.c_Ngood_min = int(cfg['Control light curve cut settings']['Ngood_min'])
			print(f'# Bound for an epoch\'s maximum chi-square: {self.c_x2_max}')
			print(f'# Bound for an epoch\'s maximum abs(flux/dflux) ratio: {self.stn_max}')
			print(f'# Bound for an epoch\'s maximum number of clipped control measurements: {self.c_Nclip_max}')
			print(f'# Bound for an epoch\'s minimum number of good control measurements: {self.c_Ngood_min}')

		self.averaging = args.average 
		if self.averaging:
			print(f'\nAveraging and cutting bad days: {self.averaging}')
			self.flux2mag_sigmalimit = int(cfg['Input/output settings']['flux2mag_sigmalimit'])
			print(f'# Sigma limit when converting flux to magnitude (magnitudes are limits when dmagnitudes are NaN): {self.flux2mag_sigmalimit}')
			self.mjd_bin_size = args.mjd_bin_size if not(args.mjd_bin_size is None) else float(cfg['Averaging settings']['mjd_bin_size'])
			print(f'# MJD bin size: {self.mjd_bin_size} days')
			self.keep_empty_bins = bool(cfg['Averaging settings']['keep_empty_bins'])
			print(f'# Keep empty bins and store as NaN in averaged light curve: {self.keep_empty_bins}')
			
			self.g_Nclip_max = int(cfg['Averaging settings']['Nclip_max'])
			self.g_Ngood_min = int(cfg['Averaging settings']['Ngood_min'])
			self.g_x2_max = float(cfg['Averaging settings']['x2_max'])
			print(f'# MJD bin bounds for not flagging as bad day: ')
			print(f'## Maximum number of clipped measurements (Nclip_max): {self.g_Nclip_max}')
			print(f'## Minimum number of good measurements (Ngood_min): {self.g_Ngood_min}')
			print(f'## Maximum chi-square (x2_max): {self.g_x2_max}')

		self.plot = args.plot
		print(f'\nPlotting: {self.plot}')

	# helper function for get_baseline_regions()
	def get_Ndays(self, SN_region_index):
		return 200 if SN_region_index == 2 else 40

	# get regions of a lc where no SN flux is present
	def get_baseline_regions(self, lc, Ndays_min):
		print('# Getting region indices around SN... ')

		baseline_ix = lc.get_baseline_ix()
		tchange1 = 58417
		tchange2 = 58882

		regions = {}
		regions['t0'] = lc.pdastro.ix_inrange(colnames=['MJD'], uplim=tchange1)
		regions['t1'] = lc.pdastro.ix_inrange(colnames=['MJD'], lowlim=tchange1, uplim=tchange2)
		regions['t2'] = lc.pdastro.ix_inrange(colnames=['MJD'], lowlim=tchange2)
		regions['b_t0'] = AandB(regions['t0'], baseline_ix)
		regions['b_t1'] = AandB(regions['t1'], baseline_ix)
		regions['b_t2'] = AandB(regions['t2'], baseline_ix)

		for region_index in range(0,3):
			if len(regions['t%d'%region_index]) > 0:
				print('## TEMPLATE REGION t%d MJD RANGE: %0.2f - %0.2f' % (region_index, lc.pdastro.t.loc[regions['t%d'%region_index][0],'MJD'], lc.pdastro.t.loc[regions['t%d'%region_index][-1],'MJD']))
			else:
				print('## TEMPLATE REGION t%d MJD RANGE: not found' % region_index)

		# find region SN starts in 
		SN_region_index = None
		if lc.discdate <= tchange1:
			SN_region_index = 0
		elif lc.discdate > tchange1 and lc.discdate <= tchange2:
			SN_region_index = 1
		elif lc.discdate > tchange2:
			SN_region_index = 2
		if SN_region_index is None:
			raise RuntimeError('## ERROR: Something went wrong--could not find region with SN discovery date!')
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
			if len(regions['b_t%d'%region_index]) > 0:
				print('## TEMPLATE REGION t%d BASELINE MJD RANGE: %0.2f - %0.2f' % (region_index, lc.pdastro.t.loc[regions['b_t%d'%region_index][0],'MJD'], lc.pdastro.t.loc[regions['b_t%d'%region_index][-1],'MJD']))
			else:
				print('## TEMPLATE REGION t%d BASELINE MJD RANGE: not found' % region_index)

		# check to make sure baseline flux is still consistent by getting median of first and last halves of affected region
		first_i = regions['b_t%d'%adjust_region_index][0]
		mid_i   = regions['b_t%d'%adjust_region_index][int(len(regions['b_t%d'%adjust_region_index])/2)]
		last_i  = regions['b_t%d'%adjust_region_index][-1]
		median1 = np.median(lc.pdastro.t.loc[lc.pdastro.ix_inrange(colnames=['MJD'], lowlim=lc.pdastro.t.loc[first_i,'MJD'], uplim=lc.pdastro.t.loc[mid_i,'MJD']), 'uJy'])
		median2 = np.median(lc.pdastro.t.loc[lc.pdastro.ix_inrange(colnames=['MJD'], lowlim=lc.pdastro.t.loc[mid_i+1,'MJD'], uplim=lc.pdastro.t.loc[last_i,'MJD']), 'uJy'])
		print(f'## Checking that baseline flux is consistent throughout adjusted region...\n## Median of first half: {median1:0.2f}\n## Median of second half: {median2:0.2f}')

		lc.corrected_baseline_ix = np.concatenate([regions['b_t0'], regions['b_t1'], regions['b_t2']])
		lc.during_sn_ix = AnotB(lc.pdastro.getindices(), lc.corrected_baseline_ix)
		
		return regions, lc

	# correct control light curves for atlas template changes at mjd=58417,58882 
	# get median of same baseline regions as SN, then apply to entire region
	def controls_correct_for_template(self, lc, control_index, regions, region_index):
		b_goodx2_i = lc.lcs[control_index].ix_inrange(colnames=['chi/N'], uplim=5)

		lowlim = lc.pdastro.t.loc[regions[f'b_t{region_index}'][0], 'MJD'] #lc.lcs[control_index].t.loc[regions[f'b_t{region_index}'][0],"MJD"]
		uplim= lc.pdastro.t.loc[regions[f'b_t{region_index}'][-1], 'MJD'] #lc.lcs[control_index].t.loc[regions[f'b_t{region_index}'][-1],"MJD"]
		b_region_i = lc.lcs[control_index].ix_inrange(colnames=['MJD'], lowlim=lowlim, uplim=uplim, exclude_uplim=True)

		if len(b_region_i) > 0:
			#print(f'#### Adjusting for template change in region b_t{region_index} from {lowlim:0.2f}-{uplim:0.2f}...')
			#print(f'#### Baseline median before: {np.median(lc.lcs[control_index].t.loc[b_region_i,"uJy"])}')
						
			if len(AandB(b_region_i,b_goodx2_i)) > 0:
				median = np.median(lc.lcs[control_index].t.loc[AandB(b_region_i,b_goodx2_i),'uJy'])
			else:
				median = np.median(lc.lcs[control_index].t.loc[b_region_i,'uJy'])

			lowlim = lc.pdastro.t.loc[regions[f't{region_index}'][0], 'MJD']
			uplim = lc.pdastro.t.loc[regions[f't{region_index}'][-1], 'MJD']
			t_region_i = lc.lcs[control_index].ix_inrange(colnames=['MJD'], lowlim=lowlim, uplim=uplim, exclude_uplim=True)

			#print(f'#### Subtracting median {median:0.1f} uJy of baseline flux with chi-square ≤ 5 from light curve flux due to potential flux in the template...')
			lc.lcs[control_index].t.loc[t_region_i,'uJy'] -= median
			#print(f'#### Baseline median now: {np.median(lc.lcs[control_index].t.loc[b_region_i,"uJy"])}')

		#else:
			#print(f'#### No baseline region for region b_t{region_index}, skipping...')

		return lc

	# correct for atlas template changes at mjd=58417,58882 
	# more info here: https://fallingstar-data.com/forcedphot/faq/
	def correct_for_template(self, lc):
		print('\nCorrecting for potential flux in template due to template changes at MJD=58417,58882...')
		
		# automatically define baseline regions according to discovery date
		regions, lc = self.get_baseline_regions(lc, Ndays_min=6)

		# get indices of measurements with x2<=5 so that when getting median, use these indices if possible
		b_goodx2_i = lc.pdastro.ix_inrange(colnames=['chi/N'], uplim=5, indices=lc.corrected_baseline_ix)

		# for each region, adjust for template change by subtracting median of that region's baseline flux
		for region_index in range(0,3):
			region_i = regions[f'b_t{region_index}']
			if len(region_i) > 0:
				print(f'# Adjusting for template change in region b_t{region_index} from {lc.pdastro.t.loc[region_i[0],"MJD"]:0.2f}-{lc.pdastro.t.loc[region_i[-1],"MJD"]:0.2f}...')
				print(f'## Baseline median before: {np.median(lc.pdastro.t.loc[region_i,"uJy"])}')
				
				if len(AandB(region_i,b_goodx2_i)) > 0:
					median = np.median(lc.pdastro.t.loc[AandB(region_i,b_goodx2_i),'uJy'])
				else:
					median = np.median(lc.pdastro.t.loc[region_i,'uJy'])

				print(f'## Subtracting median {median:0.1f} uJy of baseline flux with chi-square ≤ 5 from light curve flux due to potential flux in the template...')
				lc.pdastro.t.loc[regions['t%d'%region_index],'uJy'] -= median
				print(f'## Baseline median now: {np.median(lc.pdastro.t.loc[region_i,"uJy"])}')

				if self.controls:
					print(f'## Correcting control light curves for potential flux in template...')
					for control_index in lc.lcs: #range(1,self.num_controls+1):
						#print(f'### Control index: {control_index}')
						lc = self.controls_correct_for_template(lc, control_index, regions, region_index)
			else:
				print(f'# No baseline region for region b_t{region_index}, skipping...')

		return lc

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
			#print('Dropping extra columns: ',dropcols)
			lc.pdastro.t.drop(columns=dropcols,inplace=True)

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
	  
		b_good_i = lc.pdastro.ix_inrange(colnames=['uJy/duJy'],lowlim=-self.stn_bound,uplim=self.stn_bound,indices=indices)
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

		return data

	# get the optimal contamination and loss cuts according to percentages in limcuts table
	def get_limcuts(self, lc, limcuts):
		contam_cut = None
		loss_cut = None
		contam_case = None
		loss_case = None

		sortby_loss = limcuts.t.iloc[(limcuts.t['Ploss']).argsort()].reset_index()
		min_loss = sortby_loss.loc[0,'Ploss']
		max_loss = sortby_loss.loc[len(sortby_loss)-1,'Ploss']
		# if all loss below lim, loss_cut is min cut
		if min_loss < self.loss_lim and max_loss < self.loss_lim:
			loss_case = 'below lim'
			loss_cut = limcuts.t.loc[0,'PSF Chi-Square Cut']
		else:
			# else if all loss above lim, loss_cut is min cut with min% loss
			if min_loss > self.loss_lim and max_loss > self.loss_lim:
				loss_case = 'above lim'
				a = np.where(limcuts.t['Ploss'] == min_loss)[0]
				b = limcuts.t.iloc[a]
				c = b.iloc[(b['PSF Chi-Square Cut']).argsort()].reset_index()
				loss_cut = c.loc[0,'PSF Chi-Square Cut']
			# else if loss crosses lim at some point, loss_cut is min cut with max% loss <= loss_lim
			else:
				loss_case = 'crosses lim'
				valid_cuts = sortby_loss[sortby_loss['Ploss'] <= self.loss_lim]
				a = np.where(limcuts.t['Ploss'] == valid_cuts.loc[len(valid_cuts)-1,'Ploss'])[0]
				# sort by cuts
				b = limcuts.t.iloc[a]
				c = b.iloc[(b['PSF Chi-Square Cut']).argsort()].reset_index()
				# get midpoint of loss1 and loss2 (two points on either side of lim)
				loss1_i = np.where(limcuts.t['PSF Chi-Square Cut'] == c.loc[0,'PSF Chi-Square Cut'])[0][0]
				if limcuts.t.loc[loss1_i,'Ploss'] == self.loss_lim:
					loss_cut = limcuts.t.loc[loss1_i,'PSF Chi-Square Cut']
				else:
					loss2_i = loss1_i - 1
					x = np.array([limcuts.t.loc[loss1_i,'PSF Chi-Square Cut'], limcuts.t.loc[loss2_i,'PSF Chi-Square Cut']])
					contam_y = np.array([limcuts.t.loc[loss1_i,'Pcontamination'], limcuts.t.loc[loss2_i,'Pcontamination']])
					loss_y = np.array([limcuts.t.loc[loss1_i,'Ploss'], limcuts.t.loc[loss2_i,'Ploss']])
					contam_line = np.polyfit(x,contam_y,1)
					loss_line = np.polyfit(x,loss_y,1)
					loss_cut = (self.loss_lim-loss_line[1])/loss_line[0]

		sortby_contam = limcuts.t.iloc[(limcuts.t['Pcontamination']).argsort()].reset_index()
		min_contam = sortby_contam.loc[0,'Pcontamination']
		max_contam = sortby_contam.loc[len(sortby_contam)-1,'Pcontamination']
		# if all contam below lim, contam_cut is max cut
		if min_contam < self.contam_lim and max_contam < self.contam_lim:
			contam_case = 'below lim'
			contam_cut = limcuts.t.loc[len(limcuts.t)-1,'PSF Chi-Square Cut']
		else:
			# else if all contam above lim, contam_cut is max cut with min% contam
			if min_contam > self.contam_lim and max_contam > self.contam_lim:
				contam_case = 'above lim'
				a = np.where(limcuts.t['Pcontamination'] == min_contam)[0]
				b = limcuts.t.iloc[a]
				c = b.iloc[(b['PSF Chi-Square Cut']).argsort()].reset_index()
				contam_cut = c.loc[len(c)-1,'PSF Chi-Square Cut']
			# else if contam crosses lim at some point, contam_cut is max cut with max% contam <= contam_lim
			else:
				contam_case = 'crosses lim'
				valid_cuts = sortby_contam[sortby_contam['Pcontamination'] <= self.contam_lim]
				a = np.where(limcuts.t['Pcontamination'] == valid_cuts.loc[len(valid_cuts)-1,'Pcontamination'])[0]
				# sort by cuts
				b = limcuts.t.iloc[a]
				c = b.iloc[(b['PSF Chi-Square Cut']).argsort()].reset_index()
				# get midpoint of contam1 and contam2 (two points on either side of lim)
				contam1_i = np.where(limcuts.t['PSF Chi-Square Cut'] == c.loc[len(c)-1,'PSF Chi-Square Cut'])[0][0]
				if limcuts.t.loc[contam1_i,'Pcontamination'] == self.contam_lim:
					contam_cut = limcuts.t.loc[contam1_i,'PSF Chi-Square Cut']
				else:
					contam2_i = contam1_i + 1
					x = np.array([limcuts.t.loc[contam1_i,'PSF Chi-Square Cut'], limcuts.t.loc[contam2_i,'PSF Chi-Square Cut']])
					contam_y = np.array([limcuts.t.loc[contam1_i,'Pcontamination'], limcuts.t.loc[contam2_i,'Pcontamination']])
					loss_y = np.array([limcuts.t.loc[contam1_i,'Ploss'], limcuts.t.loc[contam2_i,'Ploss']])
					contam_line = np.polyfit(x,contam_y,1)
					loss_line = np.polyfit(x,loss_y,1)
					contam_cut = (self.contam_lim-contam_line[1])/contam_line[0]

		loss_cut_data = self.get_limcuts_data(lc, 'loss_cut', loss_cut, loss_case)
		contam_cut_data = self.get_limcuts_data(lc, 'contam_cut', contam_cut, contam_case)

		return contam_cut_data, loss_cut_data

	# choose between contamination cut and loss cut
	def get_final_chisquare_cut(self, contam_cut, loss_cut, contam_case, loss_case):
		case1 = loss_case == 'below lim' or contam_case == 'below lim'
		case2 = loss_case == 'above lim' or contam_case == 'above lim'
		case3 = loss_case == 'crosses lim' or contam_case == 'crosses lim'

		# case 1 and 1: final_cut = 3
		# case 1 and 2: take limit of case 2
		# case 1 and 3: take limit of case 3
		# case 2 and 2: print lims don't work
		# case 2 and 3: get_final_chisquare_cut
		# case 3 and 3: get_final_chisquare_cut

		final_cut = None
		if case1 and not case2 and not case3: # 1 and 1
			print(f'# Valid chi-square cut range from {loss_cut:0.2f} to {contam_cut:0.2f}! Setting to 3...')
			final_cut = self.min_cut
		elif case1: # 1
			if case2: # and 2
				if loss_case == 'above lim':
					print(f'# WARNING: contam_cut <= {contam_cut:0.2f} falls below limit {self.contam_lim:0.2f}%, but loss_cut >= {loss_cut:0.2f} falls above limit {self.loss_lim:0.2f}%! Setting to {loss_cut:0.2f}...')
					final_cut = loss_cut
				else:
					print(f'# WARNING: loss_cut <= {loss_cut:0.2f} falls below limit {self.loss_lim:0.2f}%, but contam_cut >= {contam_cut:0.2f} falls above limit {self.contam_lim:0.2f}%! Setting to {contam_cut:0.2f}...')
					final_cut = contam_cut
			else: # and 3
				if loss_case == 'crosses lim':
					print(f'# contam_cut <= {contam_cut:0.2f} falls below limit {self.contam_lim:0.2f}% and loss_cut >= {loss_cut:0.2f} crosses limit {self.loss_lim:0.2f}%, setting to {loss_cut:0.2f}...')
					final_cut = loss_cut
				else:
					print(f'# loss_cut <= {loss_cut:0.2f} falls below limit {self.loss_lim:0.2f}% and contam_cut >= {contam_cut:0.2f} crosses limit {self.contam_lim:0.2f}%, setting to {contam_cut:0.2f}...')
					final_cut = contam_cut
		elif case2 and not case3: # 2 and 2
			print(f'# ERROR: chi-square loss_cut >= {loss_cut:0.2f} and contam_cut <= {contam_cut:0.2f} both fall above respective limits {self.loss_lim:0.2f}% and {self.contam_lim:0.2f}%! Try setting less strict limits. Setting final cut to nan.')
			final_cut = np.nan
		else: # 2 and 3 or 3 and 3
			if loss_cut > contam_cut:
				print(f'# WARNING: chi-square loss_cut >= {loss_cut:0.2f} and contam_cut <= {contam_cut:0.2f} do not overlap!')
				if self.lim_to_prioritize == 'contam_lim':
					print(f'# Prioritizing {self.lim_to_prioritize} and setting to {contam_cut:0.2f}...')
					final_cut = contam_cut
				else:
					print(f'# Prioritizing {self.lim_to_prioritize} and setting to {loss_cut:0.2f}...')
					final_cut = loss_cut
			else:
				print(f'# Valid chi-square cut range from {loss_cut:0.2f} to {contam_cut:0.2f}!')
				if self.lim_to_prioritize == 'contam_lim':
					print(f'# Prioritizing {self.lim_to_prioritize} and setting to {loss_cut:0.2f}...')
					final_cut = loss_cut
				else:
					print(f'# Prioritizing {self.lim_to_prioritize} and setting to {contam_cut:0.2f}...')
					final_cut = contam_cut
		return final_cut

	# apply chi-square cut to SN light curve and update mask column with flag
	def apply_chisquare_cut(self, args, lc, plot=None):
		print('\nNow applying dynamic chi-square cut...')

		if self.chisquare_cut is None:
			limcuts = self.get_limcuts_table(lc)
			contam_cut_data, loss_cut_data = self.get_limcuts(lc, limcuts)

			if args.plot:
				if plot is None:
					raise RuntimeError('Plot object not passed to apply_chisquare_cut()!')
				plot.plot_limcuts(limcuts, contam_cut_data['cut'], loss_cut_data['cut'], self.min_cut, self.max_cut)

			print(f'# Contamination cut according to given contam_limit, with {contam_cut_data["Pcontamination"]:0.2f}% contamination and {contam_cut_data["Ploss"]:0.2f}% loss: {contam_cut_data["cut"]:0.2f}')
			if contam_cut_data['case'] == 'above lim':
				print(f'## WARNING: Contamination cut not possible with contamination <= contam_lim {self.contam_lim:0.2f}%!')
			print(f'# Loss cut according to given loss_limit, with {loss_cut_data["Pcontamination"]:0.2f}% contamination and {loss_cut_data["Ploss"]:0.2f}% loss: {loss_cut_data["cut"]:0.2f}')
			if loss_cut_data['case'] == 'above lim':
				print(f'## WARNING: Loss cut not possible with loss <= loss_lim {self.loss_lim:0.2f}%!')

			final_cut = self.get_final_chisquare_cut(contam_cut_data['cut'], loss_cut_data['cut'], contam_cut_data['case'], loss_cut_data['case'])

			if np.isnan(final_cut):
				raise RuntimeError('\n# ERROR: Final suggested chi-square cut could not be determined according to given contamination and loss limits. We suggest resetting your limits in atlas_lc_settings.ini.')
			else:
				if final_cut == contam_cut_data['cut']:
					Pcontamination = contam_cut_data['Pcontamination']
					Ploss = contam_cut_data['Ploss']
				else:
					Pcontamination = loss_cut_data['Pcontamination']
					Ploss = loss_cut_data['Ploss']
				print(f'# Final suggested chi-square cut is {final_cut:0.2f}, with {Pcontamination:0.2f}% contamination and {Ploss:0.2f}% loss.')
				if (Pcontamination > self.contam_lim):
					print(f'## WARNING: Final cut\'s contamination {Pcontamination:0.2f}% exceeds {self.contam_lim:0.2f}%!')
				if (Ploss > self.loss_lim):
					print(f'## WARNING: Final cut\'s loss {Ploss:0.2f}% exceeds loss_lim {self.loss_lim:0.2f}%!')
		else:
			final_cut = self.chisquare_cut
			print(f'Chi-square cut set to {final_cut} manually by user, overriding dynamic chi-square cut')

		# update mask column with final chi-square cut
		cut_ix = lc.pdastro.ix_inrange(colnames=['chi/N'], lowlim=final_cut, exclude_lowlim=True)
		lc.update_mask_col(self.flags['chisquare'], cut_ix)
		print(f'# Total percent of data flagged: {100*len(cut_ix)/len(lc.pdastro.getindices()):0.2f}%')

		if args.plot:
			return lc, final_cut, plot
		else:
			return lc, final_cut

	# apply chi-square cut to SN light curve and update mask column with flag
	def apply_uncertainty_cut(self, lc):
		print('\nNow applying uncertainty cut...')

		# update mask column with uncertainty cut
		cut_ix = lc.pdastro.ix_inrange(colnames=['duJy'], lowlim=self.uncertainty_cut, exclude_lowlim=True)
		lc.update_mask_col(self.flags['uncertainty'], cut_ix)
		print(f'# Total percent of data flagged: {100*len(cut_ix)/len(lc.pdastro.getindices()):0.2f}%')

		return lc

	# make sure that for every SN measurement, we have corresponding control light curve 
	# measurements at that MJD
	def verify_mjds(self, lc):
		print('# Making sure SN and control light curve MJDs match up exactly...')
		# sort SN lc by MJD
		mjd_sorted_i = lc.pdastro.ix_sort_by_cols('MJD')
		lc.pdastro.t = lc.pdastro.t.loc[mjd_sorted_i]
		sn_sorted = lc.pdastro.t.loc[mjd_sorted_i,'MJD'].to_numpy()

		for control_index in lc.lcs: #range(1,self.num_controls+1):
			# sort control light curves by MJD
			mjd_sorted_i = lc.lcs[control_index].ix_sort_by_cols('MJD')
			control_sorted = lc.lcs[control_index].t.loc[mjd_sorted_i,'MJD'].to_numpy()
			
			# compare control light curve to SN light curve and, if out of agreement, fix
			if (len(sn_sorted) != len(control_sorted)) or (np.array_equal(sn_sorted, control_sorted) is False):
				print('## MJDs out of agreement for control light curve %03d, fixing...' % control_index)

				mjds_onlysn = AnotB(sn_sorted, control_sorted)
				mjds_onlycontrol = AnotB(control_sorted, sn_sorted)

				# for the MJDs only in SN, add row with that MJD to control light curve, with all values of other columns NaN
				if len(mjds_onlysn) > 0:
					print('### Adding %d NaN rows to control light curve...' % len(mjds_onlysn))
					for mjd in mjds_onlysn:
						lc.lcs[control_index].newrow({'MJD':mjd,'Mask':0})
				
				# remove indices of rows in control light curve for which there is no MJD in the SN lc
				if len(mjds_onlycontrol) > 0:
					print('### Removing %d control light curve rows without matching SN rows...' % len(mjds_onlycontrol))
					indices2skip = []
					for mjd in mjds_onlycontrol:
						ix = lc.lcs[control_index].ix_equal('MJD',mjd)
						if len(ix)!=1:
							raise RuntimeError(f'### Couldn\'t find MJD={mjd} in column MJD, but should be there!')
						indices2skip.extend(ix)
					indices = AnotB(lc.lcs[control_index].getindices(),indices2skip)
				else:
					indices = lc.lcs[control_index].getindices()
				
				ix_sorted = lc.lcs[control_index].ix_sort_by_cols('MJD',indices=indices)
				lc.lcs[control_index].t = lc.lcs[control_index].t.loc[ix_sorted]
		
		return lc 

	def get_control_stats(self, lc):
		print('# Calculating control light curve statistics...')

		# construct arrays for control lc data
		uJy = np.full((self.num_controls, len(lc.pdastro.t['MJD'])), np.nan)
		duJy = np.full((self.num_controls, len(lc.pdastro.t['MJD'])), np.nan)
		Mask = np.full((self.num_controls, len(lc.pdastro.t['MJD'])), 0, dtype=np.int32)
		
		for control_index in lc.lcs: #range(1,self.num_controls+1):
			if (len(lc.lcs[control_index].t) != len(lc.pdastro.t['MJD'])) or (np.array_equal(lc.pdastro.t['MJD'], lc.lcs[control_index].t['MJD']) is False):
				raise RuntimeError(f'## sERROR: SN lc not equal to control lc for control_index {control_index}! Rerun or debug verify_mjds().')
			else:
				uJy[control_index-1,:] = lc.lcs[control_index].t['uJy']
				duJy[control_index-1,:] = lc.lcs[control_index].t['duJy']
				Mask[control_index-1,:] = lc.lcs[control_index].t['Mask']

		c2_param2columnmapping = lc.pdastro.intializecols4statparams(prefix='c2_',format4outvals='{:.2f}',skipparams=['converged','i'])

		for index in range(uJy.shape[-1]):
			pda4MJD = pdastrostatsclass()
			pda4MJD.t['uJy'] = uJy[1:,index]
			pda4MJD.t['duJy'] = duJy[1:,index]
			pda4MJD.t['Mask'] = np.bitwise_and(Mask[1:,index], self.flags['chisquare']|self.flags['uncertainty'])
			
			pda4MJD.calcaverage_sigmacutloop('uJy',noisecol='duJy',maskcol='Mask',maskval=(self.flags['chisquare']|self.flags['uncertainty']),verbose=1,Nsigma=3.0,median_firstiteration=True)
			lc.pdastro.statresults2table(pda4MJD.statparams, c2_param2columnmapping, destindex=index) 

		return lc

	def print_control_flag_stats(self, lc):
		print('# Control light curve cut results:')
		print('## Length of SN light curve: %d' % len(lc.pdastro.t))
		print('## Percent of data above x2_max bound: %0.2f%%' % (100*len(lc.pdastro.ix_masked('Mask',maskval=self.flags['controls_x2']))/len(lc.pdastro.t)))
		print('## Percent of data above stn_max bound: %0.2f%%' % (100*len(lc.pdastro.ix_masked('Mask',maskval=self.flags['controls_stn']))/len(lc.pdastro.t)))
		print('## Percent of data above Nclip_max bound: %0.2f%%' % (100*len(lc.pdastro.ix_masked('Mask',maskval=self.flags['controls_Nclip']))/len(lc.pdastro.t)))
		print('## Percent of data below Ngood_min bound: %0.2f%%' % (100*len(lc.pdastro.ix_masked('Mask',maskval=self.flags['controls_Ngood']))/len(lc.pdastro.t)))
		print('## Total percent of data flagged as bad: %0.2f%%' % (100*len(lc.pdastro.ix_masked('Mask',maskval=self.flags['controls_bad']))/len(lc.pdastro.t)))
		print('## Total percent of data flagged as questionable (not masked with control light curve flags but Nclip > 0): %0.2f%%' % (100*len(lc.pdastro.ix_masked('Mask',maskval=self.flags['controls_questionable']))/len(lc.pdastro.t)))

	# apply control light curve cut to SN light curve and update mask column with flag
	def apply_control_cut(self, lc):
		print('\nNow applying control light curve cut...')

		# clear any previous flags in control light curves' 'Mask' columns
		for control_index in lc.lcs: #range(1,self.num_controls+1):
			lc.lcs[control_index].t['Mask'] = 0

		lc = self.verify_mjds(lc)
		lc = self.get_control_stats(lc)

		print('# Flagging SN light curve based on control light curve statistics...')
		lc.pdastro.t['c2_abs_stn'] = lc.pdastro.t['c2_mean']/lc.pdastro.t['c2_mean_err']

		# flag measurements according to given bounds
		flag_x2_i = lc.pdastro.ix_inrange(colnames=['c2_X2norm'], lowlim=self.c_x2_max, exclude_lowlim=True)
		lc.update_mask_col(self.flags['controls_x2'], flag_x2_i)
		flag_stn_i = lc.pdastro.ix_inrange(colnames=['c2_abs_stn'], lowlim=self.stn_max, exclude_lowlim=True)
		lc.update_mask_col(self.flags['controls_stn'], flag_stn_i)
		flag_nclip_i = lc.pdastro.ix_inrange(colnames=['c2_Nclip'], lowlim=self.c_Nclip_max, exclude_lowlim=True)
		lc.update_mask_col(self.flags['controls_Nclip'], flag_nclip_i)
		flag_ngood_i = lc.pdastro.ix_inrange(colnames=['c2_Ngood'], uplim=self.c_Ngood_min, exclude_uplim=True)
		lc.update_mask_col(self.flags['controls_Ngood'], flag_ngood_i)

		# update mask column with control light curve cut on any measurements flagged according to given bounds
		zero_Nclip_i = lc.pdastro.ix_equal('c2_Nclip', 0)
		unmasked_i = lc.pdastro.ix_unmasked('Mask', maskval=self.flags['controls_x2']|self.flags['controls_stn']|self.flags['controls_Nclip']|self.flags['controls_Ngood'])
		lc.update_mask_col(self.flags['controls_questionable'], AnotB(unmasked_i,zero_Nclip_i))
		lc.update_mask_col(self.flags['controls_bad'], AnotB(lc.pdastro.getindices(),unmasked_i))
		self.print_control_flag_stats(lc)

		return lc

	# average the SN light curve
	def average_lc(self, lc, avglc):
		print(f'\nNow averaging SN light curve...')

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
				if fluxstatparams['Ngood'] < self.g_Ngood_min:
					is_bad = True
				if fluxstatparams['Nclip'] > self.g_Nclip_max:
					is_bad = True
				if not(fluxstatparams['X2norm'] is None) and fluxstatparams['X2norm'] > self.g_x2_max:
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

		print('# Total percent of data flagged: %0.2f%%' % (100*len(lc.pdastro.ix_masked('Mask',maskval=self.flags['avg_badday']))/len(avglc.pdastro.t)))

		return lc, avglc

	# output text file with information on cuts and flags
	def add_to_readme(self, f, lc, filt, final_cut):
		f.write(f'\n\n## FILTER: {filt}')

		if self.chisquares:
			f.write(f'\n\n### Chi-square cut')
			f.write(f'\nWe flag measurements with a chi-square (column name "chi/N") value above {final_cut:0.2f} with hex value {hex(self.flags["chisquare"])}.')

		if self.uncertainties:
			f.write(f'\n\n### Uncertainty cut')
			f.write(f'\nWe flag measurements with an uncertainty (column name "duJy") value above {self.uncertainty_cut:0.2f} with hex value {hex(self.flags["uncertainty"])}.')
		
		if self.controls:
			f.write(f'\n\n### Control light curve cut')
			f.write(f'\nThe control light curve cut examines each SN epoch and its corresponding control light curve measurements at that epoch, applies a 3-sigma-clipped average, calculates statistics, and then cuts bad epochs based on those returned statistics.')
			f.write(f'\n\nFor the given epoch, we flag the SN measurement for which the returned control statistics fulfill any of the following criteria with the hex value {hex(self.flags["controls_bad"])}: ')
			f.write(f'\n\t- A returned chi-square > {self.c_x2_max:0.2f}')
			f.write(f'\n\t- A returned abs(flux/dflux) > {self.stn_max:0.2f}')
			f.write(f'\n\t- Number of clipped/"bad" measurements in the 3σ-clipped average > {self.c_Nclip_max}')
			f.write(f'\n\t- Number of used/"good" measurements in the 3σ-clipped average < {self.c_Ngood_min}')

		f.write(f'\n\nAfter any or all of these three cuts are applied, the light curves are saved with the new Mask column.')

		if self.averaging:
			f.write(f'\n\n### Averaging light curves and cutting bad days')
			f.write(f'\nFor each MJD bin of size {self.mjd_bin_size:0.2f} day(s), we calculate the 3σ-clipped average of any SN measurements falling within that bin and use that average as our flux for that bin. However, out of all exposures within this MJD bin, only measurements not cut in the previous methods are averaged in the 3σ-clipped average cut. (The exception to this statement would be the case that all 4 measurements are cut in previous methods; in this case, they are averaged anyway and flagged as a bad bin.')
			f.write(f'\n\nThen we flag any measurements in the SN light curve for the given epoch for which statistics fulfill any of the following criteria with the hex value {hex(self.flags["avg_badday"])}: ') 
			f.write(f'\n\t- A returned chi-square > {self.g_x2_max}')
			f.write(f'\n\t- Number of measurements averaged < {self.g_Ngood_min}')
			f.write(f'\n\t- Number of measurements clipped > {self.g_Nclip_max}')
			f.write(f'\n\nThe averaged light curves are then saved in a new file with the MJD bin size added to the filename.')

		return f

	def get_lc_data(self, lc, snlist_index):
		if snlist_index == -1:
			lc.get_tns_data(self.tns_api_key)

			# add row to self.snlist
			self.snlist.newrow({'tnsname':lc.tnsname, 
								'ra':lc.ra, 
								'dec':lc.dec, 
								'discovery_date':lc.discdate, 
								'closebright_ra':np.nan, 
								'closebright_dec':np.nan})
			
			snlist_index = len(self.snlist.t) - 1
		
		lc.ra = self.snlist.t.loc[snlist_index,'ra']
		lc.dec = self.snlist.t.loc[snlist_index,'dec']
		lc.discdate = self.snlist.t.loc[snlist_index,'discovery_date']
		print(f'# RA: {lc.ra}, Dec: {lc.dec}, discovery date: {lc.discdate}')

		return lc, snlist_index

	def cut_loop(self):
		args = self.define_args().parse_args()
		self.load_settings(args)

		for obj_index in range(0,len(args.tnsnames)):
			print(f'\nCOMMENCING LOOP FOR SN {args.tnsnames[obj_index]}')

			f = open(f'{self.output_dir}/{args.tnsnames[obj_index]}/README.md','w+')
			f.write(f"# SN {args.tnsnames[obj_index]} Light Curve Cleaning and Averaging")

			if args.plot:
				plot = plot_atlas_lc(tnsname=args.tnsnames[obj_index], output_dir=self.output_dir, args=args, flags=self.flags)

			snlist_index = -1
			snlist_ix = self.snlist.ix_equal(colnames=['tnsname'],val=args.tnsnames[obj_index])
			# check if SN information exists in snlist.txt
			if len(snlist_ix) > 0:
				if len(snlist_ix > 1):
					# drop duplicate rows
					self.snlist.t.drop(snlist_ix[1:])
				snlist_index = snlist_ix[0]

			for filt in ['o','c']:
				print(f'\nFILTER SET: {filt}')
				lc = atlas_lc(tnsname=args.tnsnames[obj_index])
				lc.load(filt, self.output_dir, num_controls=self.num_controls)
				lc, snlist_index = self.get_lc_data(lc, snlist_index) #lc.get_tns_data(self.tns_api_key) # TO DO: NO NEED TO REPEAT THIS FOR EACH FILTER--ADD SNLIST.TXT?

				lc = self.drop_extra_columns(lc)
				lc = self.correct_for_template(lc)

				# clear 'Mask' column
				lc.pdastro.t['Mask'] = 0

				# add flux/dflux column
				print('Adding uJy/duJy column to light curve...')
				lc.pdastro.t['uJy/duJy'] = lc.pdastro.t['uJy']/lc.pdastro.t['duJy']
				lc.pdastro.t = lc.pdastro.t.replace([np.inf, -np.inf], np.nan)

				if args.plot:
					plot.set(lc=lc, filt=filt)
					plot.plot_og_lc()

				if self.chisquares:
					if args.plot:
						lc, final_cut, plot = self.apply_chisquare_cut(args, lc, plot)
						plot.plot_chisquare_cut()
					else:
						lc, final_cut = self.apply_chisquare_cut(args, lc)

				if self.uncertainties:
					lc = self.apply_uncertainty_cut(lc)
					if args.plot:
						plot.plot_uncertainty_cut()

				if self.controls:
					lc = self.apply_control_cut(lc)
					if args.plot:
						plot.plot_controls_cut()

				if args.plot and (self.chisquares or self.uncertainties or self.controls):
					plot.plot_all_cuts()

				if self.averaging:
					avglc = atlas_lc(tnsname=lc.tnsname, is_averaged=True, mjd_bin_size=self.mjd_bin_size)
					avglc.discdate = lc.discdate

					lc, avglc = self.average_lc(lc, avglc)

					# save averaged light curve
					avglc.save(self.output_dir, filt=filt, overwrite=self.overwrite)

					# plot averaged light curve
					if args.plot:
						plot.set(lc=avglc, filt=filt)
						plot.plot_averaged_lc()

				# drop extra control lc cut columns and save lc with new 'Mask' column
				lc = self.drop_extra_columns(lc)
				lc.save(self.output_dir, filt=filt, overwrite=self.overwrite)

				f = self.add_to_readme(f, lc, filt, final_cut)

			f.close()

			if args.plot:
				plot.save()

		# save snlist.txt with any new rows
		print(f'Saving SN list at {self.snlist_filename}')
		self.snlist.write(self.snlist_filename)

if __name__ == "__main__":
	clean_atlas_lc = clean_atlas_lc()
	clean_atlas_lc.cut_loop()
