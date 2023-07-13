#!/usr/bin/env python
"""
@author: Sofia Rest
"""

import sys, argparse, configparser, re, os
from copy import deepcopy
import pandas as pd
import numpy as np

from pdastro import pdastrostatsclass, AandB, AnotB
from atlas_lc import atlas_lc
from plot_atlas_lc import plot_atlas_lc
from asym_gaussian import gauss2lc

class clean_atlas_lc():
	def __init__(self):
		# credentials
		self.tns_api_key = None
		self.tns_id = None
		self.bot_name = None

		# input/output
		self.output_dir = None
		self.snlist_filename = None
		self.snlist = None
		self.overwrite = True
		self.num_controls = 0

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

		# uncertainty cut
		self.uncertainties = False 
		self.uncertainty_cut = None

		# estimating true uncertainties
		self.estimate_true_uncertainties = False
		self.estimate_true_uncertainties_chisquare_cut = None

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

		# control light curve cut
		self.controls = False
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

		# detecting pre-SN bumps
		self.detect_bumps = False
		self.gaussian_sigma = None
		self.appmags = None
		self.start_mjd = None
		self.end_mjd = None
	
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

		parser.add_argument('-b', '--detect_bumps', default=False, action='store_true', help='apply rolling gaussian weighted sum to flux/dflux in order to amplify possible precursor bumps')
		parser.add_argument('--sim_gaussian', nargs=3, default=None, help=('comma-separated peakMJD list, peak_appmag, gaussian_sigma: add a gaussian at peakMJD with a peak apparent magnitude of peak_appmag and a sigma of gaussian_sigma in days'))

		parser.add_argument('-p', '--plot', default=False, action='store_true', help='plot each cut and save into PDF file')
		parser.add_argument('--xlim_lower', type=float, default=None, help='if plotting, manually set lower x axis limit to a certain MJD')
		parser.add_argument('--xlim_upper', type=float, default=None, help='if plotting, manually set upper x axis limit to a certain MJD')
		parser.add_argument('--ylim_lower', type=float, default=None, help='if plotting, manually set lower y axis limit to a certain uJy')
		parser.add_argument('--ylim_upper', type=float, default=None, help='if plotting, manually set upper y axis limit to a certain uJy')

		parser.add_argument('-f','--cfg_filename', default='params.ini', type=str, help='file name of ini file with settings for this class')
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
		self.tns_id = cfg['TNS credentials']['tns_id']
		self.bot_name = cfg['TNS credentials']['bot_name']
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
		self.num_controls = int(cfg['Control light curve settings']['num_controls'])
		print(f'Number of control light curves: {self.num_controls}')
		self.plot = args.plot
		print(f'Plotting: {self.plot}')

		self.estimate_true_uncertainties = cfg['True uncertainties estimation settings']['estimate_true_uncertainties']=='True'
		self.estimate_true_uncertainties_chisquare_cut = float(cfg['True uncertainties estimation settings']['estimate_true_uncertainties_chisquare_cut'])
		if self.estimate_true_uncertainties:
			print(f'Estimating true uncertainties set to {self.estimate_true_uncertainties} with preliminary chi-square cut at {self.estimate_true_uncertainties_chisquare_cut:0.2f}')

		self.chisquares = args.chisquares
		if self.chisquares:
			print(f'\nChi-square cut: {self.chisquares}')
			try:
				self.chisquare_cut = float(cfg['Chi-square cut settings']['override_cut'])
				print(f'# Overriding dynamic chi-square cut with manual cut of x2 = {self.chisquare_cut}')
			except:
				self.stn_bound = float(cfg['Chi-square cut settings']['stn_bound'])
				self.min_cut = int(cfg['Chi-square cut settings']['min_cut'])
				self.max_cut = int(cfg['Chi-square cut settings']['max_cut'])
				self.cut_step = int(cfg['Chi-square cut settings']['cut_step'])
				self.contam_lim = float(cfg['Chi-square cut settings']['contamination_limit'])
				self.loss_lim = float(cfg['Chi-square cut settings']['loss_limit'])
				self.lim_to_prioritize = cfg['Chi-square cut settings']['limit_to_prioritize']
				if not(self.lim_to_prioritize == 'loss' or self.lim_to_prioritize == 'contamination'):
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
			self.keep_empty_bins = cfg['Averaging settings']['keep_empty_bins']=='True'
			print(f'# Keep empty bins and store as NaN in averaged light curve: {self.keep_empty_bins}')
			
			self.g_Nclip_max = int(cfg['Averaging settings']['Nclip_max'])
			self.g_Ngood_min = int(cfg['Averaging settings']['Ngood_min'])
			self.g_x2_max = float(cfg['Averaging settings']['x2_max'])
			print(f'# MJD bin bounds for not flagging as bad day: ')
			print(f'## Maximum number of clipped measurements (Nclip_max): {self.g_Nclip_max}')
			print(f'## Minimum number of good measurements (Ngood_min): {self.g_Ngood_min}')
			print(f'## Maximum chi-square (x2_max): {self.g_x2_max}')

		self.detect_bumps = args.detect_bumps
		if self.detect_bumps:
			print(f'\nDetecting pre-SN bumps: {self.detect_bumps}')
			self.apply_to_controls = cfg['Detecting bumps settings']['apply_to_controls']=='True'
			print(f'# Applying to control light curves in order to establish detection limit: {self.apply_to_controls}')
			self.gaussian_sigma = float(cfg['Detecting bumps settings']['gaussian_sigma'])
			print(f'# Searching for pre-SN bumps with a sigma of {self.gaussian_sigma:0.2f} days')
			
			# simulated gaussian settings
			if not(args.sim_gaussian is None):
				print(f'# Adding simulated gaussian pre-SN bump to SN light curve: True')
				if ',' in args.sim_gaussian[1]:
					self.appmags = args.sim_gaussian[1].split(',')
					print(f'## Multiple magnitudes input: {self.appmags}')
				else:
					self.appmags = [args.sim_gaussian[1]]
					print(f'## Only one magnitude input: {self.appmags}')
			
			# custom MJD range for bump detection
			if not(cfg['Detecting bumps settings']['start_mjd'] == 'None'):
				self.start_mjd = float(cfg['Detecting bumps settings']['start_mjd'])
				print('# Will start bump detection at MJD={self.start_mjd}')
			if not(cfg['Detecting bumps settings']['end_mjd'] == 'None'):
				print('# Will end bump detection at MJD={self.end_mjd}')
				self.end_mjd = float(cfg['Detecting bumps settings']['end_mjd'])

	# helper function for get_baseline_regions()
	def get_Ndays(self, SN_region_index):
		return 200 if SN_region_index == 2 else 40

	def get_SNstart_region(self, discdate, tchange1, tchange2):
		# find region SN starts in 
		SN_region_index = None
		if discdate <= tchange1:
			SN_region_index = 0
		elif discdate > tchange1 and discdate <= tchange2:
			SN_region_index = 1
		elif discdate > tchange2:
			SN_region_index = 2
		if SN_region_index is None:
			raise RuntimeError('## ERROR: Something went wrong--could not find region with SN discovery date!')
		else:
			print('## SN discovery date located in template region t%d' % SN_region_index)
		return SN_region_index


	# get regions of a lc where no SN flux is present
	def get_baseline_regions(self, lc, Ndays_min):
		print('# Getting region indices around SN... ')

		baseline_ix = lc.get_baseline_ix()
		tchange1 = 58417
		tchange2 = 58882

		regions = {}
		regions['t0'] = lc.lcs[0].ix_inrange(colnames=['MJD'], uplim=tchange1)
		regions['t1'] = lc.lcs[0].ix_inrange(colnames=['MJD'], lowlim=tchange1, uplim=tchange2)
		regions['t2'] = lc.lcs[0].ix_inrange(colnames=['MJD'], lowlim=tchange2)
		regions['b_t0'] = AandB(regions['t0'], baseline_ix)
		regions['b_t1'] = AandB(regions['t1'], baseline_ix)
		regions['b_t2'] = AandB(regions['t2'], baseline_ix)

		found_region_ix = []
		for region_index in range(0,3):
			if len(regions['t%d'%region_index]) > 0:
				print('## TEMPLATE REGION t%d MJD RANGE: %0.2f - %0.2f' % (region_index, lc.lcs[0].t.loc[regions['t%d'%region_index][0],'MJD'], lc.lcs[0].t.loc[regions['t%d'%region_index][-1],'MJD']))
				found_region_ix.append(region_index)
			else:
				print('## TEMPLATE REGION t%d MJD RANGE: not found' % region_index)
		
		# cannot do flux correction?
		if len(found_region_ix) <= 1:
			print('WARNING: At least 2 template regions do not contain any data. Therefore, flux could not be corrected according to the ATLAS reference template changes. Skipping...')
			
			# try to get baseline flux
			SN_region_index = self.get_SNstart_region(lc.discdate, tchange1, tchange2) # find region with SN in it
			if SN_region_index in found_region_ix: # SN discovery date is in one of the found regions
				# get last Ndays days for baseline
				lc.corrected_baseline_ix = lc.lcs[0].ix_inrange(colnames=['MJD'],
																lowlim=lc.lcs[0].t.loc[regions['t%d'%SN_region_index][-1],'MJD'] - self.get_Ndays(SN_region_index),
																uplim=lc.lcs[0].t.loc[regions['t%d'%SN_region_index][-1],'MJD'])
			elif len(found_region_ix) == 1: # one template region has data
				# get all of found template for baseline
				lc.corrected_baseline_ix = regions[f't{found_region_ix[0]}']
			else: # no template region has data, so basically no data... this is pretty bad but hopefully will not come to this
				lc.corrected_baseline_ix = baseline_ix # will be None
			lc.during_sn_ix = AnotB(lc.lcs[0].getindices(), lc.corrected_baseline_ix)

			return None, lc

		SN_region_index = self.get_SNstart_region(lc.discdate, tchange1, tchange2)

		# for region with tail end of the SN, get last Ndays days and classify as baseline
		adjust_region_index = SN_region_index
		if adjust_region_index < 2 and len(regions['b_t%d'%adjust_region_index]) >= Ndays_min:
			adjust_region_index += 1
		if len(regions['b_t%d'%adjust_region_index]) < Ndays_min:
			print('## Getting baseline flux for template region t%d by obtaining last %d days of region... ' % (adjust_region_index, self.get_Ndays(adjust_region_index)))
			regions['b_t%d'%adjust_region_index] = lc.lcs[0].ix_inrange(colnames=['MJD'],
																		lowlim=lc.lcs[0].t.loc[regions['t%d'%adjust_region_index][-1],'MJD'] - self.get_Ndays(adjust_region_index),
																		uplim=lc.lcs[0].t.loc[regions['t%d'%adjust_region_index][-1],'MJD'])
		if adjust_region_index < 1: regions['b_t1'] = regions['t1']
		if adjust_region_index < 2: regions['b_t2'] = regions['t2']

		for region_index in range(0,3):
			if len(regions['b_t%d'%region_index]) > 0:
				print('## TEMPLATE REGION t%d BASELINE MJD RANGE: %0.2f - %0.2f' % (region_index, lc.lcs[0].t.loc[regions['b_t%d'%region_index][0],'MJD'], lc.lcs[0].t.loc[regions['b_t%d'%region_index][-1],'MJD']))
			else:
				print('## TEMPLATE REGION t%d BASELINE MJD RANGE: not found' % region_index)

		# check to make sure baseline flux is still consistent by getting median of first and last halves of affected region
		first_i = regions['b_t%d'%adjust_region_index][0]
		mid_i   = regions['b_t%d'%adjust_region_index][int(len(regions['b_t%d'%adjust_region_index])/2)]
		last_i  = regions['b_t%d'%adjust_region_index][-1]
		median1 = np.median(lc.lcs[0].t.loc[lc.lcs[0].ix_inrange(colnames=['MJD'], lowlim=lc.lcs[0].t.loc[first_i,'MJD'], uplim=lc.lcs[0].t.loc[mid_i,'MJD']), 'uJy'])
		median2 = np.median(lc.lcs[0].t.loc[lc.lcs[0].ix_inrange(colnames=['MJD'], lowlim=lc.lcs[0].t.loc[mid_i+1,'MJD'], uplim=lc.lcs[0].t.loc[last_i,'MJD']), 'uJy'])
		print(f'## Checking that baseline flux is consistent throughout adjusted region...\n## Median of first half: {median1:0.2f}\n## Median of second half: {median2:0.2f}')

		lc.corrected_baseline_ix = np.concatenate([regions['b_t0'], regions['b_t1'], regions['b_t2']])
		lc.during_sn_ix = AnotB(lc.lcs[0].getindices(), lc.corrected_baseline_ix)
		
		return regions, lc

	# correct control light curves for atlas template changes at mjd=58417,58882 
	def controls_correct_for_template(self, lc, control_index, regions, region_index):
		goodx2_i = lc.lcs[control_index].ix_inrange(colnames=['chi/N'], uplim=5)

		# get indices of control lc that match up with SN's baseline region
		#lowlim = lc.lcs[0].t.loc[regions[f'b_t{region_index}'][0], 'MJD'] 
		#uplim = lc.lcs[0].t.loc[regions[f'b_t{region_index}'][-1], 'MJD']
		
		# get indices of target template region
		tchange1 = 58417
		tchange2 = 58882
		lowlim = None
		uplim = None
		if region_index == 0:
			uplim = tchange1
		elif region_index == 1:
			lowlim = tchange1
			uplim = tchange2
		elif region_index == 2:
			lowlim = tchange2
		region_i = lc.lcs[control_index].ix_inrange(colnames=['MJD'], lowlim=lowlim, uplim=uplim, exclude_uplim=True)

		if len(region_i) > 0:
			# get median of template region
			if len(AandB(region_i,goodx2_i)) > 0:
				median = np.median(lc.lcs[control_index].t.loc[AandB(region_i,goodx2_i),'uJy'])
			else:
				median = np.median(lc.lcs[control_index].t.loc[region_i,'uJy'])

			#lowlim = lc.lcs[0].t.loc[regions[f't{region_index}'][0], 'MJD']
			#uplim = lc.lcs[0].t.loc[regions[f't{region_index}'][-1], 'MJD']
			#t_region_i = lc.lcs[control_index].ix_inrange(colnames=['MJD'], lowlim=lowlim, uplim=uplim, exclude_uplim=True)
			#lc.lcs[control_index].t.loc[t_region_i,'uJy'] -= median

			# subtract median from entire template region of control lc
			lc.lcs[control_index].t.loc[region_i,'uJy'] -= median

		return lc

	# correct for atlas template changes at mjd=58417,58882 
	# more info here: https://fallingstar-data.com/forcedphot/faq/
	def correct_for_template(self, lc):
		print('\nCorrecting for potential flux in template due to template changes at MJD=58417,58882...')
		output = ''
		
		# automatically define baseline regions according to discovery date
		regions, lc = self.get_baseline_regions(lc, Ndays_min=6)
		if regions is None:
			# TODO: ADD SECTION TO README
			output += 'Not enough data found in at least 2 template regions. Could not correct for template changes.'
			return lc, output

		# get indices of measurements with x2<=5 so that when getting median, use these indices if possible
		b_goodx2_i = lc.lcs[0].ix_inrange(colnames=['chi/N'], uplim=5, indices=lc.corrected_baseline_ix)

		# for each region, adjust for template change by subtracting median of that region's baseline flux
		for region_index in range(0,3):
			region_i = regions[f'b_t{region_index}']
			if len(region_i) > 0:
				print(f'# Adjusting for template change in region b_t{region_index} from {lc.lcs[0].t.loc[region_i[0],"MJD"]:0.2f}-{lc.lcs[0].t.loc[region_i[-1],"MJD"]:0.2f}...')
				print(f'## Baseline median before: {np.median(lc.lcs[0].t.loc[region_i,"uJy"])}')
				
				if len(AandB(region_i,b_goodx2_i)) > 0:
					median = np.median(lc.lcs[0].t.loc[AandB(region_i,b_goodx2_i),'uJy'])
				else:
					median = np.median(lc.lcs[0].t.loc[region_i,'uJy'])

				print(f'## Subtracting median {median:0.1f} uJy of baseline flux with chi-square ≤ 5 from light curve flux due to potential flux in the template...')
				lc.lcs[0].t.loc[regions['t%d'%region_index],'uJy'] -= median
				print(f'## Baseline median now: {np.median(lc.lcs[0].t.loc[region_i,"uJy"])}')
				output += f'\nCorrection applied to baseline region {region_index}: {median:0.1f} uJy subtracted'

				# control lc correction
				if self.controls:
					print(f'## Correcting control light curves for potential flux in template...')
					for control_index in range(1,self.num_controls+1):
						lc = self.controls_correct_for_template(lc, control_index, regions, region_index)
			else:
				print(f'# No baseline region for region b_t{region_index}, skipping...')

		return lc, output

	# drop mask column and any added columns from previous iterations
	def drop_extra_columns(self, lc, control_index=0):
		dropcols=[]

		for col in ['Noffsetlc', 'uJy/duJy', '__tmp_SN', 'SNR', 'SNRsum', 'SNRsumnorm', 'SNRsim', 'SNRsimsum']:
			if col in lc.lcs[control_index].t.columns:
				dropcols.append(col)
		for col in lc.lcs[control_index].t.columns:
			if re.search('^c\d_',col): 
				dropcols.append(col)

		# drop any extra columns
		if len(dropcols)>0: 
			#print('Dropping extra columns: ',dropcols)
			lc.lcs[control_index].t.drop(columns=dropcols,inplace=True)

		return lc

	def apply_true_uncertainties(self, lc):
		output = ''
		for control_index in range(self.num_controls+1):
			if control_index == 0:
				print('\nNow estimating true uncertainties for SN light curve...')

				if len(lc.corrected_baseline_ix) <= 0:
					print('WARNING: No available baseline flux! Cannot proceed with true uncertainties estimation. Skipping...')
					output += '\nNot enough available baseline flux! Could not proceed with true uncertainties estimation.'
					for control_index in range(self.num_controls+1):
						lc.dflux_colnames[control_index] = 'duJy'
					return lc, output, False

				clean_ix = AandB(lc.lcs[0].ix_unmasked('Mask',maskval=self.flags['uncertainty']), lc.lcs[0].ix_inrange(colnames=['chi/N'],uplim=self.estimate_true_uncertainties_chisquare_cut,exclude_uplim=True))
				clean_ix = AandB(lc.corrected_baseline_ix, clean_ix)
			else:
				print(f'Now estimating true uncertainties for control light curve {control_index:03d}...')
				clean_ix = AandB(lc.lcs[control_index].ix_unmasked('Mask',maskval=self.flags['uncertainty']), lc.lcs[control_index].ix_inrange(colnames=['chi/N'],uplim=self.estimate_true_uncertainties_chisquare_cut,exclude_uplim=True))

			lc.lcs[control_index].calcaverage_sigmacutloop('uJy', indices=clean_ix, Nsigma=3.0, median_firstiteration=True, verbose=1)
			#if control_index == 0: print(lc.lcs[control_index].statparams)
			sigma_true_typical = lc.lcs[control_index].statparams['stdev']

			median_dflux = np.median(lc.lcs[control_index].t.loc[clean_ix, 'duJy'])

			if sigma_true_typical > median_dflux:
				print(f'# True typical uncertainty {sigma_true_typical:0.2f} greater the current median uncertainty {median_dflux:0.2f}. Proceeding with true uncertainties estimation...')

				# for following cuts, use updated uncertainty column for this light curve
				lc.dflux_colnames[control_index] = 'duJy_new'

				# add extra noise source to current noise in new column
				sigma_extra = np.sqrt(sigma_true_typical*sigma_true_typical - median_dflux)
				print(f'# Sigma extra calculated: {sigma_extra:0.4f}')
				#print('# Adding sigma_extra {sigma_extra:0.4f} to new duJy column...')
				lc.lcs[control_index].t['duJy_new'] = np.nan
				lc.lcs[control_index].t['duJy_new'] = np.sqrt(lc.lcs[control_index].t['duJy']**2 + sigma_extra**2)

				#lc.lcs[control_index].t['uJy/duJy_new'] = lc.lcs[control_index].t['uJy']/lc.lcs[control_index].t['duJy_new']

				if control_index == 0:
					output += f' An extra noise of sigma {sigma_extra:0.4f} was added to the uncertainties of the SN light curve and copied to the "duJy_new" column.'
			
				return lc, output, True
			else:
				print(f'# True typical uncertainty less than or equal to the current median uncertainty. Skipping true uncertainties estimation.')
				return lc, output, False


	# apply uncertainty cut to SN light curve and update mask column with flag
	def apply_uncertainty_cut(self, lc):
		print('\nNow applying uncertainty cut...')

		# update SN mask column with final chi-square cut and apply same cut to control light curves
		for control_index in range(self.num_controls+1):
			cut_ix = lc.lcs[control_index].ix_inrange(colnames=['duJy'], lowlim=self.uncertainty_cut, exclude_lowlim=True)
			lc.update_mask_col(self.flags['uncertainty'], cut_ix, control_index=control_index)
			if control_index == 0:
				s = f'Total percent of data flagged: {100*len(cut_ix)/len(lc.lcs[0].getindices()):0.2f}%'
				output = f'\n\n{s}.'
				print(f'# {s}')

		return lc, output

	# for a range of chi-square cuts, determine contamination, loss, and other pecentages
	def get_limcuts_table(self, lc, indices=None):
		limcuts = pdastrostatsclass(columns=['PSF Chi-Square Cut', 'N', 'Ngood', 'Nbad', 'Nkept', 'Ncut', 'Ngood,kept', 'Ngood,cut', 'Nbad,kept', 'Nbad,cut',
											  'Pgood,kept', 'Pgood,cut', 'Pbad,kept', 'Pbad,cut', 'Ngood,kept/Ngood', 'Ploss', 'Pcontamination',
											  'Nbad,cut 3<stn<=5', 'Nbad,cut 5<stn<=10', 'Nbad,cut 10<stn', 'Nbad,kept 3<stn<=5', 'Nbad,kept 5<stn<=10', 'Nbad,kept 10<stn'])

		if indices is None:
			indices = lc.corrected_baseline_ix
		
		# define good and bad baseline measurement indices according to abs(uJy/duJy) bound
		b_good_i = lc.lcs[0].ix_inrange(colnames=['uJy/duJy'], lowlim=-self.stn_bound, uplim=self.stn_bound, indices=indices)
		b_bad_i = AnotB(indices, b_good_i)

		# for different x2 cuts decreasing from 50
		for cut in range(self.min_cut, self.max_cut+1, self.cut_step):
			# define kept baseline measurement indices according to current chi-square cut
			b_kept_i = lc.lcs[0].ix_inrange(colnames=['chi/N'], uplim=cut, indices=indices)
			b_cut_i = AnotB(indices, b_kept_i)

			if 100*(len(b_kept_i)/len(indices)) < 10:
 				# less than 10% of measurements kept, so this chi-square cut not valid
				print(f'# At cut {cut}, less than 10% of measurements are kept ({100*(len(b_kept_i)/len(indices)):0.2f}% kept)--skipping...')
				continue
			else: 
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
								len(AandB(AandB(b_bad_i,b_cut_i), lc.lcs[0].ix_inrange(colnames=['uJy/duJy'], lowlim=-3, uplim=5, exclude_lowlim=True))), # Nbad,cut 3<stn<=5
								len(AandB(AandB(b_bad_i,b_cut_i), lc.lcs[0].ix_inrange(colnames=['uJy/duJy'], lowlim=5, uplim=10, exclude_lowlim=True))), # Nbad,cut 5<stn<=10
								len(AandB(AandB(b_bad_i,b_cut_i), lc.lcs[0].ix_inrange(colnames=['uJy/duJy'], lowlim=10, exclude_lowlim=True))), # Nbad,cut 10<stn 
								len(AandB(AandB(b_bad_i,b_kept_i), lc.lcs[0].ix_inrange(colnames=['uJy/duJy'], lowlim=-3, uplim=5, exclude_lowlim=True))), # Nbad,kept 3<stn<=5
								len(AandB(AandB(b_bad_i,b_kept_i), lc.lcs[0].ix_inrange(colnames=['uJy/duJy'], lowlim=5, uplim=10, exclude_lowlim=True))), # Nbad,kept 5<stn<=10
								len(AandB(AandB(b_bad_i,b_kept_i), lc.lcs[0].ix_inrange(colnames=['uJy/duJy'], lowlim=10, exclude_lowlim=True))), # Nbad,kept 10<stn 
								]], columns=['PSF Chi-Square Cut', 'N', 'Ngood', 'Nbad', 'Nkept', 'Ncut', 'Ngood,kept', 'Ngood,cut', 'Nbad,kept', 'Nbad,cut',
											 'Pgood,kept', 'Pgood,cut', 'Pbad,kept', 'Pbad,cut', 'Ngood,kept/Ngood', 'Ploss', 'Pcontamination',
											 'Nbad,cut 3<stn<=5', 'Nbad,cut 5<stn<=10', 'Nbad,cut 10<stn', 'Nbad,kept 3<stn<=5', 'Nbad,kept 5<stn<=10', 'Nbad,kept 10<stn'])
			limcuts.t = pd.concat([limcuts.t,df],ignore_index=True)
		return limcuts 

	# get contamination and loss plus other percentages and information for a certain chi-square cut 
	def get_limcuts_data(self, lc, name, cut, case):
		indices = lc.corrected_baseline_ix
	  
		b_good_i = lc.lcs[0].ix_inrange(colnames=['uJy/duJy'],lowlim=-self.stn_bound,uplim=self.stn_bound,indices=indices)
		b_bad_i = AnotB(indices, b_good_i)
		b_kept_i = lc.lcs[0].ix_inrange(colnames=['chi/N'],uplim=cut,indices=indices)
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
					contam_cut = (self.contam_lim - contam_line[1])/contam_line[0]

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
			print(f'# Valid chi-square cut range from {loss_cut:0.2f} to {contam_cut:0.2f}! Setting to {self.min_cut:0.2f}...')
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
				plot.plot_limcuts(limcuts, contam_cut_data['cut'], loss_cut_data['cut'], self.contam_lim, self.loss_lim, self.min_cut, self.max_cut)

			print(f'# Contamination cut according to given contam_limit, with {contam_cut_data["Pcontamination"]:0.2f}% contamination and {contam_cut_data["Ploss"]:0.2f}% loss: {contam_cut_data["cut"]:0.2f}')
			if contam_cut_data['case'] == 'above lim':
				print(f'## WARNING: Contamination cut not possible with contamination <= contam_lim {self.contam_lim:0.2f}%!')
			print(f'# Loss cut according to given loss_limit, with {loss_cut_data["Pcontamination"]:0.2f}% contamination and {loss_cut_data["Ploss"]:0.2f}% loss: {loss_cut_data["cut"]:0.2f}')
			if loss_cut_data['case'] == 'above lim':
				print(f'## WARNING: Loss cut not possible with loss <= loss_lim {self.loss_lim:0.2f}%!')

			final_cut = self.get_final_chisquare_cut(contam_cut_data['cut'], loss_cut_data['cut'], contam_cut_data['case'], loss_cut_data['case'])

			if np.isnan(final_cut):
				raise RuntimeError('\n# ERROR: Final suggested chi-square cut could not be determined according to given contamination and loss limits. We suggest resetting your limits in params.ini.')
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
		
			output = f'\n\t- The cut optimized according to the given contamination limit of {self.contam_lim:0.2f}% was {contam_cut_data["cut"]:0.2f}, with a contamination of {contam_cut_data["Pcontamination"]:0.2f}% and a loss of {contam_cut_data["Ploss"]:0.2f}%.'
			output += f'\n\t- The cut optimized according to the given loss limit of {self.loss_lim:0.2f}% was {loss_cut_data["cut"]:0.2f}, with a contamination of {loss_cut_data["Pcontamination"]:0.2f}% and a loss of {loss_cut_data["Ploss"]:0.2f}%.'
		else:
			final_cut = self.chisquare_cut

			print(f'Chi-square cut set to {final_cut} manually by user, overriding dynamic chi-square cut')
			output = f'Cut was manually set to {final_cut} by user, overriding dynamic chi-square cut.'

		# update SN mask column with final chi-square cut and apply same cut to control light curves
		for control_index in range(self.num_controls+1):
			cut_ix = lc.lcs[control_index].ix_inrange(colnames=['chi/N'], lowlim=final_cut, exclude_lowlim=True)
			lc.update_mask_col(self.flags['chisquare'], cut_ix, control_index=control_index)
			if control_index == 0:
				s = f'Total percent of data flagged: {100*len(cut_ix)/len(lc.lcs[0].getindices()):0.2f}%'
				output += f'\n\n{s}.'
				print(f'# {s}')

		if args.plot:
			return lc, final_cut, output, plot
		else:
			return lc, final_cut, output

	# make sure that for every SN measurement, we have corresponding control light curve 
	# measurements at that MJD
	def verify_mjds(self, lc):
		print('# Making sure SN and control light curve MJDs match up exactly...')
		# sort SN lc by MJD
		mjd_sorted_i = lc.lcs[0].ix_sort_by_cols('MJD')
		lc.lcs[0].t = lc.lcs[0].t.loc[mjd_sorted_i]
		sn_sorted = lc.lcs[0].t.loc[mjd_sorted_i,'MJD'].to_numpy()

		for control_index in range(1,self.num_controls+1):
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
					print('### Removing %d control light curve row(s) without matching SN row(s)...' % len(mjds_onlycontrol))
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
		uJy = np.full((self.num_controls, len(lc.lcs[0].t['MJD'])), np.nan)
		duJy = np.full((self.num_controls, len(lc.lcs[0].t['MJD'])), np.nan)
		Mask = np.full((self.num_controls, len(lc.lcs[0].t['MJD'])), 0, dtype=np.int32)
		
		for control_index in range(1,self.num_controls+1):
			if (len(lc.lcs[control_index].t) != len(lc.lcs[0].t['MJD'])) or (np.array_equal(lc.lcs[0].t['MJD'], lc.lcs[control_index].t['MJD']) is False):
				raise RuntimeError(f'## sERROR: SN lc not equal to control lc for control_index {control_index}! Rerun or debug verify_mjds().')
			else:
				uJy[control_index-1,:] = lc.lcs[control_index].t['uJy']
				duJy[control_index-1,:] = lc.lcs[control_index].t[lc.dflux_colnames[control_index]]
				Mask[control_index-1,:] = lc.lcs[control_index].t['Mask']

		c2_param2columnmapping = lc.lcs[0].intializecols4statparams(prefix='c2_',format4outvals='{:.2f}',skipparams=['converged','i'])

		for index in range(uJy.shape[-1]):
			pda4MJD = pdastrostatsclass()
			pda4MJD.t['uJy'] = uJy[0:,index]
			pda4MJD.t['duJy'] = duJy[0:,index]
			pda4MJD.t['Mask'] = np.bitwise_and(Mask[0:,index], self.flags['chisquare']|self.flags['uncertainty'])
			
			pda4MJD.calcaverage_sigmacutloop('uJy',noisecol='duJy',maskcol='Mask',maskval=(self.flags['chisquare']|self.flags['uncertainty']),verbose=1,Nsigma=3.0,median_firstiteration=True)
			lc.lcs[0].statresults2table(pda4MJD.statparams, c2_param2columnmapping, destindex=index) 

		return lc

	def print_control_flag_stats(self, lc):
		print('# Control light curve cut results:')
		print('## Length of SN light curve: %d' % len(lc.lcs[0].t))
		print('## Percent of data above x2_max bound: %0.2f%%' % (100*len(lc.lcs[0].ix_masked('Mask',maskval=self.flags['controls_x2']))/len(lc.lcs[0].t)))
		print('## Percent of data above stn_max bound: %0.2f%%' % (100*len(lc.lcs[0].ix_masked('Mask',maskval=self.flags['controls_stn']))/len(lc.lcs[0].t)))
		print('## Percent of data above Nclip_max bound: %0.2f%%' % (100*len(lc.lcs[0].ix_masked('Mask',maskval=self.flags['controls_Nclip']))/len(lc.lcs[0].t)))
		print('## Percent of data below Ngood_min bound: %0.2f%%' % (100*len(lc.lcs[0].ix_masked('Mask',maskval=self.flags['controls_Ngood']))/len(lc.lcs[0].t)))
		
		s = 'Total percent of data flagged as bad: %0.2f%%' % (100*len(lc.lcs[0].ix_masked('Mask',maskval=self.flags['controls_bad']))/len(lc.lcs[0].t))
		print(f'## {s}')
		output = f'\n\n {s}.'
		s = 'Total percent of data flagged as questionable (not masked with control light curve flags but Nclip > 0): %0.2f%%' % (100*len(lc.lcs[0].ix_masked('Mask',maskval=self.flags['controls_questionable']))/len(lc.lcs[0].t))
		print(f'## {s}')
		output += f'\n {s}.'

		return output

	# apply control light curve cut to SN light curve and update mask column with flag
	def apply_control_cut(self, lc):
		print('\nNow applying control light curve cut...')

		# clear any previous flags in control light curves' 'Mask' columns
		#for control_index in range(1,self.num_controls+1):
			#lc.lcs[control_index].t['Mask'] = 0

		#lc = self.verify_mjds(lc)
		lc = self.get_control_stats(lc)

		print('# Flagging SN light curve based on control light curve statistics...')
		lc.lcs[0].t['c2_abs_stn'] = lc.lcs[0].t['c2_mean']/lc.lcs[0].t['c2_mean_err']

		# flag measurements according to given bounds
		flag_x2_i = lc.lcs[0].ix_inrange(colnames=['c2_X2norm'], lowlim=self.c_x2_max, exclude_lowlim=True)
		lc.update_mask_col(self.flags['controls_x2'], flag_x2_i)
		flag_stn_i = lc.lcs[0].ix_inrange(colnames=['c2_abs_stn'], lowlim=self.stn_max, exclude_lowlim=True)
		lc.update_mask_col(self.flags['controls_stn'], flag_stn_i)
		flag_nclip_i = lc.lcs[0].ix_inrange(colnames=['c2_Nclip'], lowlim=self.c_Nclip_max, exclude_lowlim=True)
		lc.update_mask_col(self.flags['controls_Nclip'], flag_nclip_i)
		flag_ngood_i = lc.lcs[0].ix_inrange(colnames=['c2_Ngood'], uplim=self.c_Ngood_min, exclude_uplim=True)
		lc.update_mask_col(self.flags['controls_Ngood'], flag_ngood_i)

		# update mask column with control light curve cut on any measurements flagged according to given bounds
		zero_Nclip_i = lc.lcs[0].ix_equal('c2_Nclip', 0)
		unmasked_i = lc.lcs[0].ix_unmasked('Mask', maskval=self.flags['controls_x2']|self.flags['controls_stn']|self.flags['controls_Nclip']|self.flags['controls_Ngood'])
		lc.update_mask_col(self.flags['controls_questionable'], AnotB(unmasked_i,zero_Nclip_i))
		lc.update_mask_col(self.flags['controls_bad'], AnotB(lc.lcs[0].getindices(),unmasked_i))

		# copy over SN's control cut flags to control light curve 'Mask' column
		flags_arr = np.full(lc.lcs[0].t['Mask'].shape, (self.flags['controls_bad']|
														self.flags['controls_questionable']|
														self.flags['controls_x2']|
														self.flags['controls_stn']|
														self.flags['controls_Nclip']|
														self.flags['controls_Ngood']))
		flags_to_copy = np.bitwise_and(lc.lcs[0].t['Mask'], flags_arr)
		for control_index in range(1,self.num_controls+1):
			lc.lcs[control_index].t['Mask'] = lc.lcs[control_index].t['Mask'].astype(np.int32)
			if len(lc.lcs[control_index].t) < 1:
				continue
			elif len(lc.lcs[control_index].t) == 1:
				lc.lcs[control_index].t.loc[0,'Mask']= int(lc.lcs[control_index].t.loc[0,'Mask']) | flags_to_copy
			else:
				lc.lcs[control_index].t['Mask'] = np.bitwise_or(lc.lcs[control_index].t['Mask'], flags_to_copy)
		
		output = self.print_control_flag_stats(lc)

		return lc, output

	def average_lc(self, lc, avglc, control_index=0):
		if control_index == 0:
			print(f'\nNow averaging SN light curve...')
		else:
			print(f'Now averaging control light curve {control_index:03d}...')

		avglc.lcs[control_index] = pdastrostatsclass()

		mjd = int(np.amin(lc.lcs[control_index].t['MJD']))
		mjd_max = int(np.amax(lc.lcs[control_index].t['MJD']))+1

		good_i = lc.lcs[control_index].ix_unmasked('Mask', maskval=self.flags['chisquare']|self.flags['uncertainty']|self.flags['controls_bad'])

		while mjd <= mjd_max:
			range_i = lc.lcs[control_index].ix_inrange(colnames=['MJD'], lowlim=mjd, uplim=mjd+self.mjd_bin_size, exclude_uplim=True)
			range_good_i = AandB(range_i,good_i)

			# add new row to averaged light curve if keep_empty_bins or any measurements present
			#if self.keep_empty_bins or len(range_i) >= 1:
			new_row = {'MJDbin':mjd+0.5*self.mjd_bin_size, 'Nclip':0, 'Ngood':0, 'Nexcluded':len(range_i)-len(range_good_i), 'Mask':0}
			avglc_index = avglc.lcs[control_index].newrow(new_row)
			
			# if no measurements present, flag or skip over day
			if len(range_i) < 1:
				#if self.keep_empty_bins:
				avglc.update_mask_col(self.flags['avg_badday'], [avglc_index], control_index=control_index)
				mjd += self.mjd_bin_size
				continue
			
			# if no good measurements, average values anyway and flag
			if len(range_good_i) < 1:
				# average flux
				lc.lcs[control_index].calcaverage_sigmacutloop('uJy', noisecol=lc.dflux_colnames[control_index], indices=range_i, Nsigma=3.0, median_firstiteration=True)
				fluxstatparams = deepcopy(lc.lcs[control_index].statparams)

				# get average mjd
				# TO DO: SHOULD NOISECOL HERE BE DUJY OR NONE??
				lc.lcs[control_index].calcaverage_sigmacutloop('MJD', indices=fluxstatparams['ix_good'], Nsigma=0, median_firstiteration=False)
				avg_mjd = lc.lcs[control_index].statparams['mean']

				# add row and flag
				avglc.lcs[control_index].add2row(avglc_index, {'MJD':avg_mjd, 
															   'uJy':fluxstatparams['mean'], 
															   'duJy':fluxstatparams['mean_err'], 
															   'stdev':fluxstatparams['stdev'],
															   'x2':fluxstatparams['X2norm'],
															   'Nclip':fluxstatparams['Nclip'],
															   'Ngood':fluxstatparams['Ngood'],
															   'Mask':0})
				lc.update_mask_col(self.flags['avg_badday'], range_i, control_index=control_index)
				avglc.update_mask_col(self.flags['avg_badday'], [avglc_index], control_index=control_index)

				mjd += self.mjd_bin_size
				continue
			
			# average good measurements
			lc.lcs[control_index].calcaverage_sigmacutloop('uJy', noisecol=lc.dflux_colnames[control_index], indices=range_good_i, Nsigma=3.0, median_firstiteration=True)
			fluxstatparams = deepcopy(lc.lcs[control_index].statparams)

			if fluxstatparams['mean'] is None or len(fluxstatparams['ix_good']) < 1:
				lc.update_mask_col(self.flags['avg_badday'], range_i, control_index=control_index)
				avglc.update_mask_col(self.flags['avg_badday'], [avglc_index], control_index=control_index)
				mjd += self.mjd_bin_size
				continue

			# get average mjd
			# TO DO: SHOULD NOISECOL HERE BE DUJY OR NONE??
			lc.lcs[control_index].calcaverage_sigmacutloop('MJD', noisecol=lc.dflux_colnames[control_index], indices=fluxstatparams['ix_good'], Nsigma=0, median_firstiteration=False)
			avg_mjd = lc.lcs[control_index].statparams['mean']

			# add row to averaged light curve
			avglc.lcs[control_index].add2row(avglc_index, {'MJD':avg_mjd, 
														   'uJy':fluxstatparams['mean'], 
														   'duJy':fluxstatparams['mean_err'], 
														   'stdev':fluxstatparams['stdev'],
														   'x2':fluxstatparams['X2norm'],
														   'Nclip':fluxstatparams['Nclip'],
														   'Ngood':fluxstatparams['Ngood'],
														   'Mask':0})
			
			# flag clipped measurements in lc
			if len(fluxstatparams['ix_clip']) > 0:
				lc.update_mask_col(self.flags['avg_ixclip'], fluxstatparams['ix_clip'], control_index=control_index)
			
			# if small number within this bin, flag measurements
			if len(range_good_i) < 3:
				lc.update_mask_col(self.flags['avg_smallnum'], range_good_i, control_index=control_index) # TO DO: CHANGE TO RANGE_I??
				avglc.update_mask_col(self.flags['avg_smallnum'], [avglc_index], control_index=control_index)
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
					lc.update_mask_col(self.flags['avg_badday'], range_i, control_index=control_index)
					avglc.update_mask_col(self.flags['avg_badday'], [avglc_index], control_index=control_index)

			mjd += self.mjd_bin_size
		
		# convert flux to magnitude and dflux to dmagnitude
		avglc.lcs[control_index].flux2mag('uJy','duJy','m','dm', zpt=23.9, upperlim_Nsigma=self.flux2mag_sigmalimit)

		avglc = self.drop_extra_columns(avglc)

		for col in ['Nclip','Ngood','Nexcluded','Mask']: 
			avglc.lcs[control_index].t[col] = avglc.lcs[control_index].t[col].astype(np.int32)

		return lc, avglc

	# average the SN light curve and, if necessary, control light curves
	def average(self, lc, avglc):
		for control_index in range(self.num_controls+1):
			# only average control light curves if detecting pre-SN bumps
			if (not(self.detect_bumps) or (self.detect_bumps and not(self.apply_to_controls))) and control_index > 0:
				break

			lc, avglc = self.average_lc(lc, avglc, control_index)

		s = 'Total percent of data flagged: %0.2f%%' % (100*len(lc.lcs[0].ix_masked('Mask',maskval=self.flags['avg_badday']))/len(avglc.lcs[0].t))
		print(f'# {s}')
		output = f'\n\n{s}.'

		return lc, avglc, output

	# add simulated bump if necessary and apply rolling gaussian weighted sum to light curve
	def apply_gaussian(self, avglc, control_index=0, simparams=None, filt=None):
		if self.start_mjd is None:
			self.start_mjd = avglc.lcs[control_index].t['MJDbin'].iloc[0]
		if self.end_mjd is None:
			self.end_mjd = avglc.lcs[control_index].t['MJDbin'].iloc[-1]
		ix = avglc.lcs[control_index].ix_inrange(colnames=['MJDbin'],lowlim=self.start_mjd, uplim=self.end_mjd)
		if len(ix) <= 0:
			raise RuntimeError('ERROR: specified MJD range less than or equal to 0!')
		print(f'# Applying detection to specified MJD bin range: {avglc.lcs[control_index].t["MJDbin"].iloc[ix[0]]} to  {avglc.lcs[control_index].t["MJDbin"].iloc[ix[-1]]}')

		good_ix = avglc.lcs[control_index].ix_unmasked('Mask',self.flags['avg_badday'])

		# make sure there are no lingering simulations
		dropcols=[]
		for col in ['uJysim','SNRsim','simLC','SNRsimsum']:
			if col in avglc.lcs[control_index].t.columns:
				dropcols.append(col)
		if len(dropcols) > 0:
			avglc.lcs[control_index].t.drop(columns=dropcols,inplace=True)

		avglc.lcs[control_index].t.loc[ix, 'SNR'] = 0.0
		avglc.lcs[control_index].t.loc[good_ix,'SNR'] = avglc.lcs[control_index].t.loc[good_ix,'uJy']/avglc.lcs[control_index].t.loc[good_ix,'duJy']

		if not(simparams is None):
			peakMJDs = simparams['sim_peakMJD'].split(',')
			
			# get the simulated gaussian
			mjds = avglc.lcs[control_index].t.loc[good_ix,'MJD']
			avglc.lcs[control_index].t.loc[good_ix,'uJysim'] = avglc.lcs[control_index].t.loc[good_ix,'uJy']
			avglc.lcs[control_index].t.loc[ix,'simLC'] = 0.0
			for peakMJD in peakMJDs:
				peakMJD = float(peakMJD)
				print(f'## Adding simulated gaussian at peak MJD {peakMJD:0.2f} with apparent magnitude {simparams["sim_appmag"]:0.2f}, sigma- of {simparams["sim_sigma_minus"]:0.2f}, and sigma+ of {simparams["sim_sigma_plus"]:0.2f}')

				# get simulated gaussian flux and add to light curve flux
				simflux = gauss2lc(mjds, peakMJD, simparams['sim_sigma_minus'], simparams['sim_sigma_plus'], app_mag=simparams['sim_appmag'])
				avglc.lcs[control_index].t.loc[good_ix,'uJysim'] += simflux

				# get the simulated lc for all MJDs
				simflux_all = gauss2lc(avglc.lcs[control_index].t.loc[ix,'MJDbin'], peakMJD, simparams['sim_sigma_minus'], simparams['sim_sigma_plus'], app_mag=simparams['sim_appmag'])
				avglc.lcs[control_index].t.loc[ix,'simLC'] += simflux_all

			# make sure all bad rows have SNRsim = 0.0 so they have no impact on the rolling SNRsum
			avglc.lcs[control_index].t.loc[ix,'SNRsim'] = 0.0
			# include simflux in the SNR
			avglc.lcs[control_index].t.loc[good_ix,'SNRsim'] = avglc.lcs[control_index].t.loc[good_ix,'uJysim']/avglc.lcs[control_index].t.loc[good_ix,'duJy']

		gaussian_sigma = round(self.gaussian_sigma/self.mjd_bin_size)
		windowsize = int(6*gaussian_sigma)
		halfwindowsize = int(windowsize*0.5)+1
		print(f'## Sigma (days): {self.gaussian_sigma:0.2f}; MJD bin size (days): {self.mjd_bin_size:0.2f}; sigma (bins): {gaussian_sigma:0.2f}; window size (bins): {windowsize}')

		# calculate the rolling SNR sum
		
		dataindices = np.array(range(len(avglc.lcs[control_index].t.loc[ix])) + np.full(len(avglc.lcs[control_index].t.loc[ix]), halfwindowsize))
		
		temp = pd.Series(np.zeros(len(avglc.lcs[control_index].t.loc[ix]) + 2*halfwindowsize), name='SNR', dtype=np.float64)
		temp[dataindices] = avglc.lcs[control_index].t.loc[ix,'SNR']

		SNRsum = temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=gaussian_sigma)
		avglc.lcs[control_index].t.loc[ix,'SNRsum'] = list(SNRsum[dataindices])
		
		norm_temp = pd.Series(np.zeros(len(avglc.lcs[control_index].t.loc[ix]) + 2*halfwindowsize), name='norm', dtype=np.float64)
		norm_temp[np.array(range(len(avglc.lcs[control_index].t.loc[ix])) + np.full(len(avglc.lcs[control_index].t.loc[ix]), halfwindowsize))] = np.ones(len(avglc.lcs[control_index].t.loc[ix]))
		norm_temp_sum = norm_temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=gaussian_sigma)
		avglc.lcs[control_index].t.loc[ix,'SNRsumnorm'] = list(SNRsum.loc[dataindices] / norm_temp_sum.loc[dataindices] * max(norm_temp_sum.loc[dataindices]))

		# calculate the rolling SNR sum for SNR with simflux
		if not(simparams is None):
			temp = pd.Series(np.zeros(len(avglc.lcs[control_index].t.loc[ix]) + 2*halfwindowsize), name='SNRsim', dtype=np.float64)
			temp[dataindices] = avglc.lcs[control_index].t.loc[ix,'SNRsim']
			SNRsimsum = temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=gaussian_sigma)
			avglc.lcs[control_index].t.loc[ix,'SNRsimsum'] = list(SNRsimsum.loc[dataindices])

		return avglc

	# begin output text file with information on cuts and flags
	def begin_readme(self, args, obj_index):
		f = open(f'{self.output_dir}/{args.tnsnames[obj_index]}/README.md','w+')
		f.write(f"# SN {args.tnsnames[obj_index]} Light Curve Cleaning and Averaging")
		f.write(f'\n\nThe ATLAS SN light curves are separated by filter (orange and cyan) and labelled as such in the file name. Averaged light curves contain an additional number in the file name that represents the MJD bin size used. Control light curves are located in the "controls" subdirectory and follow the same naming scheme, only with their control index added after the SN name.')
		
		f.write(f'\n\nThe following details the file names for each of the light curve versions:')
		f.write(f'\n\t- SN light curves: {args.tnsnames[obj_index]}.o.lc.txt and {args.tnsnames[obj_index]}.c.lc.txt')
		if self.averaging:
			f.write(f'\n\t- Averaged light curves: {args.tnsnames[obj_index]}.o.{self.mjd_bin_size:0.2f}days.lc.txt and {args.tnsnames[obj_index]}.c.{self.mjd_bin_size:0.2f}days.lc.txt')
		if self.controls:
			f.write(f'\n\t- Control light curves, where X=001,...,{self.num_controls:03d}: {args.tnsnames[obj_index]}_iX.o.lc.txt and {args.tnsnames[obj_index]}_iX.c.lc.txt')

		f.write(f'\n\nThe following summarizes the hex values in the "Mask" column of each light curve for each cut applied (see below sections for more information on each cut): ')
		if self.uncertainties:
			f.write(f'\n\t- Uncertainty cut: {hex(self.flags["uncertainty"])}')
		if self.chisquares:
			f.write(f'\n\t- Chi-square cut: {hex(self.flags["chisquare"])}')
		if self.controls:
			f.write(f'\n\t- Control light curve cut: {hex(self.flags["controls_bad"])}')
		if self.averaging:
			f.write(f'\n\t- Bad day (for averaged light curves): {hex(self.flags["avg_badday"])}')

		return f

	# add information about each cut to output text file
	def add_to_readme(self, f, lc, filt, final_cut=None, estimate_true_uncertainties_output=None, chisquare_output=None, uncertainty_output=None, templates_output=None, control_output=None, averaging_output=None):
		f.write(f'\n\n## FILTER: {filt}')


		f.write('\n\n### Correction for ATLAS reference template changes')
		f.write(f'\nWe take into account ATLAS\'s periodic replacement of the difference image reference templates, which may cause step discontinuities in flux. Two template changes have been recorded at MJDs 58417 and 58882.')
		f.write(templates_output)

		if self.uncertainties:
			f.write(f'\n\n### Uncertainty cut')
			f.write(f'\nWe flag measurements with an uncertainty (column name "duJy") value above {self.uncertainty_cut:0.2f} with hex value {hex(self.flags["uncertainty"])}.')
			f.write(uncertainty_output)

		if self.estimate_true_uncertainties:
			f.write(f'\n\n### Estimating true uncertainties')
			f.write(f'\nThis procedure attempts to account for an extra noise source in the data by estimating the true typical uncertainty, deriving the additional systematic uncertainty, and lastly applying this extra noise to a new uncertainty column "duJy_new". This new uncertainty column will be used in the cuts following this portion.')
			f.write(estimate_true_uncertainties_output)

		if self.chisquares:
			f.write(f'\n\n### Chi-square cut')
			f.write(f'\nWe flag measurements with a chi-square (column name "chi/N") value above {final_cut:0.2f} with hex value {hex(self.flags["chisquare"])}.')
			f.write(chisquare_output)

		if self.controls:
			f.write(f'\n\n### Control light curve cut')
			f.write(f'\nThe control light curve cut examines each SN epoch and its corresponding control light curve measurements at that epoch, applies a 3-sigma-clipped average, calculates statistics, and then cuts bad epochs based on those returned statistics.')
			f.write(f'\n\nFor the given epoch, we flag the SN measurement for which the returned control statistics fulfill any of the following criteria with the hex value {hex(self.flags["controls_bad"])}: ')
			f.write(f'\n\t- A returned chi-square > {self.c_x2_max:0.2f}')
			f.write(f'\n\t- A returned abs(flux/dflux) > {self.stn_max:0.2f}')
			f.write(f'\n\t- Number of clipped/"bad" measurements in the 3σ-clipped average > {self.c_Nclip_max}')
			f.write(f'\n\t- Number of used/"good" measurements in the 3σ-clipped average < {self.c_Ngood_min}')
			f.write(control_output)

		f.write(f'\n\n### After the uncertainty, chi-square, and control light curve cuts are applied, the light curves are resaved with the new "Mask" column.')

		if self.averaging:
			f.write(f'\n\n### Averaging light curves and cutting bad days')
			f.write(f'\nFor each MJD bin of size {self.mjd_bin_size:0.2f} day(s), we calculate the 3σ-clipped average of any SN measurements falling within that bin and use that average as our flux for that bin. However, out of all exposures within this MJD bin, only measurements not cut in the previous methods are averaged in the 3σ-clipped average cut. (The exception to this statement would be the case that all 4 measurements are cut in previous methods; in this case, they are averaged anyway and flagged as a bad bin.')
			f.write(f'\n\nThen we flag any measurements in the SN light curve for the given epoch for which statistics fulfill any of the following criteria with the hex value {hex(self.flags["avg_badday"])}: ') 
			f.write(f'\n\t- A returned chi-square > {self.g_x2_max}')
			f.write(f'\n\t- Number of measurements averaged < {self.g_Ngood_min}')
			f.write(f'\n\t- Number of measurements clipped > {self.g_Nclip_max}')
			f.write(f'\n\nThe averaged light curves are then saved in a new file with the MJD bin size added to the filename.')
			f.write(averaging_output)

		return f

	def get_lc_data(self, lc, snlist_index):
		if snlist_index == -1:
			lc.get_tns_data(self.tns_api_key, self.tns_id, self.bot_name)

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
		print(f'RA: {lc.ra}, Dec: {lc.dec}, discovery date: {lc.discdate}')

		return lc, snlist_index

	def cut_loop(self):
		args = self.define_args().parse_args()
		self.load_settings(args)

		for obj_index in range(0,len(args.tnsnames)):
			print(f'\nCOMMENCING LOOP FOR SN {args.tnsnames[obj_index]}')

			f = self.begin_readme(args, obj_index)
			
			if args.plot:
				plot = plot_atlas_lc(tnsname=args.tnsnames[obj_index], 
									 output_dir=self.output_dir, 
									 args=args, 
									 flags=self.flags)
			if self.detect_bumps:
				bumps_plot = plot_atlas_lc(tnsname=args.tnsnames[obj_index], 
										   output_dir=self.output_dir, 
										   args=args, 
										   add2filename='detect_bumps', 
										   flags=self.flags)

			# check if SN information exists in snlist.txt
			snlist_index = -1
			snlist_ix = self.snlist.ix_equal(colnames=['tnsname'],val=args.tnsnames[obj_index])
			if len(snlist_ix) > 0:
				if len(snlist_ix > 1):
					# drop duplicate rows
					self.snlist.t.drop(snlist_ix[1:])
				snlist_index = snlist_ix[0]

			for filt in ['o','c']:
				print(f'\nFILTER SET: {filt}')
				lc = atlas_lc(tnsname=args.tnsnames[obj_index])
				lc.load(self.output_dir, filt, num_controls=self.num_controls)

				lc = self.verify_mjds(lc)
				if len(lc.lcs[0].t) < 1:
					print('WARNING: Empty light curve--skipping any cuts/averaging/other...')
					continue

				lc, snlist_index = self.get_lc_data(lc, snlist_index)
				lc = self.drop_extra_columns(lc)
				lc, templates_output = self.correct_for_template(lc)
				if args.plot:
					plot.set(lc=lc, filt=filt)
					plot.plot_og_lc()

				# uncertainty cut
				uncertainty_output = None
				if self.uncertainties:
					lc, uncertainty_output = self.apply_uncertainty_cut(lc)
					if args.plot:
						plot.plot_uncertainty_cut(add2title=f'at {self.uncertainty_cut:0.2f}')

				# estimate true uncertainties
				estimate_true_uncertainties_output = None
				if self.estimate_true_uncertainties:
					lc, estimate_true_uncertainties_output, do_plot = self.apply_true_uncertainties(lc)
					if args.plot and do_plot:
						plot.plot_uncertainty_estimations()

				# add flux/dflux column
				print('Adding uJy/duJy column to light curve...')
				lc.lcs[0].t['uJy/duJy'] = lc.lcs[0].t['uJy']/lc.lcs[0].t[lc.dflux_colnames[0]]
				lc.lcs[0].t = lc.lcs[0].t.replace([np.inf, -np.inf], np.nan)

				# chi-square cut
				chisquare_output = None
				final_cut = None
				if self.chisquares:
					if args.plot:
						lc, final_cut, chisquare_output, plot = self.apply_chisquare_cut(args, lc, plot)
						plot.plot_chisquare_cut(add2title=f'at {final_cut:0.2f}')
					else:
						lc, final_cut, chisquare_output = self.apply_chisquare_cut(args, lc)

				# control light curve cut
				control_output = None
				if self.controls:
					lc, control_output = self.apply_control_cut(lc)
					if args.plot:
						plot.plot_controls_cut(self.num_controls)

				if args.plot and (self.chisquares or self.uncertainties or self.controls):
					plot.plot_all_cuts()

				# average and cut bad days
				averaging_output = None
				if self.averaging:
					avglc = atlas_lc(tnsname=lc.tnsname, is_averaged=True, mjd_bin_size=self.mjd_bin_size)
					avglc.discdate = lc.discdate

					lc, avglc, averaging_output = self.average(lc, avglc)

					# plot averaged light curve
					if args.plot:
						plot.set(lc=avglc, filt=filt)
						plot.plot_averaged_lc()

				# detect pre-SN bumps 
				if self.detect_bumps:
					print('\nNow detecting pre-SN bumps...')

					if not self.averaging:
						raise RuntimeError('ERROR: Cannot detect pre-SN bumps without averaging! Please add -g to your command in order to average light curves.')
					if self.apply_to_controls and self.num_controls <= 0:
						raise RuntimeError('ERROR: Cannot apply to control light curves without at least one control light curve! Check the num_controls field in config file.')

					if not(self.appmags is None):
						for appmag in self.appmags:
							simparams = {'sim_peakMJD':args.sim_gaussian[0],'sim_appmag':float(appmag),'sim_sigma_minus':float(args.sim_gaussian[2]),'sim_sigma_plus':float(args.sim_gaussian[2])}
							print(f'# Simulation apparent magnitude: {simparams["sim_appmag"]:0.2f} mag')
							print(f'# Simulation peak MJD(s): {simparams["sim_peakMJD"].split(",")}')
							print(f'# Simulation gaussian sigma: {simparams["sim_sigma_plus"]:0.2f} days')

							for control_index in range(self.num_controls+1):
								# only add simulated gaussian(s) to SN light curve when applying rolling sum
								if control_index == 0:
									print(f'# Applying gaussian weighted rolling sum to SN light curve...')
									avglc = self.apply_gaussian(avglc, control_index=control_index, simparams=simparams)
								else:
									print(f'# Applying gaussian weighted rolling sum to control light curve {control_index:03d}...')
									avglc = self.apply_gaussian(avglc, control_index=control_index)

							bumps_plot = plot_atlas_lc(tnsname=lc.tnsname, output_dir=self.output_dir, args=args, add2filename=f'detect_bumps.{filt}.appmag{simparams["sim_appmag"]:0.2f}', flags=self.flags)
							bumps_plot.set(lc=avglc, filt=filt)
							bumps_plot.plot_sim_bumps(simparams=simparams)
							bumps_plot.plot_snr(simparams=simparams)
							if self.apply_to_controls:
								bumps_plot.plot_all_snr(simparams=simparams)
							bumps_plot.save()
					else:
						for control_index in range(self.num_controls+1):
							if control_index == 0:
								print(f'# Applying gaussian weighted rolling sum to SN light curve...')
							else:
								print(f'# Applying gaussian weighted rolling sum to control light curve {control_index:03d}...')
							avglc = self.apply_gaussian(avglc, control_index=control_index, filt=filt)

						bumps_plot.set(lc=avglc, filt=filt)
						bumps_plot.plot_sim_bumps()
						bumps_plot.plot_snr()
						if self.apply_to_controls:
							bumps_plot.plot_all_snr()

				# drop extra control lc cut columns
				lc = self.drop_extra_columns(lc)
				# save lc with new 'Mask' column
				lc.save(self.output_dir, filt=filt, overwrite=self.overwrite)


				# drop extra columns in averaged lc
				if self.averaging:
					avglc = self.drop_extra_columns(avglc) 
					# save averaged light curve
					avglc.save(self.output_dir, filt=filt, overwrite=self.overwrite, keep_empty_bins=self.keep_empty_bins)

				f = self.add_to_readme(f, lc, filt, final_cut=final_cut, 
													uncertainty_output=uncertainty_output,
													estimate_true_uncertainties_output=estimate_true_uncertainties_output,
													templates_output=templates_output,
													chisquare_output=chisquare_output, 
													control_output=control_output,
													averaging_output=averaging_output)

			f.close()

			if args.plot:
				plot.save()
			if self.detect_bumps and self.appmags is None:
				bumps_plot.save()

		# save snlist.txt with any new rows
		filename = f'{self.output_dir}/{self.snlist_filename}'
		print(f'Saving SN list at {filename}')
		self.snlist.write(filename)

if __name__ == "__main__":
	clean_atlas_lc = clean_atlas_lc()
	clean_atlas_lc.cut_loop()
