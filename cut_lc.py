#!/usr/bin/env python
"""
Author: Sofia Rest
"""

import sys
from copy import deepcopys
import pandas as pd
import numpy as np

import sigmacut
from pdastro import pdastrostatsclass
from light_curve import light_curve

class cut_lc():
    def __init__(self):
    	self.filt = None

    	# credentials
    	self.tns_api_key = None

		# input/output
		self.input_dir = None
		self.output_dir = None
		self.overwrite = True

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

	def cut_loop(self):
		args = self.define_args().parse_args()
		self.load_settings(args)

		for obj_index in range(0,len(args.tnsnames)):
			print(f'\nCommencing cut loop for SN {lc.tnsname}')
			for filt in ['o','c']:
				self.filt = filt
				print(f'Filter set: {self.filt}')
				lc = light_curve(tnsname=args.tnsnames[obj_index])
				lc.load(self.filt, self.input_dir, num_controls=self.num_controls)
				lc.get_tns_data(self.tns_api_key)
				lc = self.correct_for_template(lc)



if __name__ == "__main__":
	cut_lc = cut_lc()
	cut_lc.cut_loop()
