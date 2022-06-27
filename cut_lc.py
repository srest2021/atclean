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
		# output
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
		self.x2_max = None
		self.stn_max = None
		self.Nclip_max = None
		self.Ngood_min = None

        
    # define command line arguments
	def define_args(self, parser=None, usage=None, conflict_handler='resolve'):
		if parser is None:
			parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)
		
		parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')

		parser.add_argument('-f','--cfg_filename', default='atlaslc.ini', type=str, help='file name of ini file with settings for this class')
		parser.add_argument('--dont_overwrite', default=False, action='store_true', help='don\'t overwrite existing file with same file name')
		
		parser.add_argument('-x', '--chisquares', default=False, action='store_true', help='apply chi-square cut')
		parser.add_argument('-u', '--uncertainties', default=False, action='store_true', help='apply uncertainty cut')
		parser.add_argument('-c', '--controls', default=False, action='store_true', help='apply control light curve cut')
		
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
			self.x2_max = cfg['Control light curve settings']['x2_max']
			self.stn_max = cfg['Control light curve settings']['stn_max']
			self.Nclip_max = cfg['Control light curve settings']['Nclip_max']
			self.Ngood_min = cfg['Control light curve settings']['Ngood_min']
			print(f'# Bound for an epoch\'s maximum chi-square: {self.x2_max}')
			print(f'# Bound for an epoch\'s maximum abs(flux/dflux) ratio: {self.stn_max}')
			print(f'# Bound for an epoch\'s maximum number of clipped control measurements: {self.Nclip_max}')
			print(f'# Bound for an epoch\'s minimum number of good control measurements: {self.Ngood_min}')

	# correct for atlas template changes at mjd=58417,58882 
	# more info here: https://fallingstar-data.com/forcedphot/faq/
	def correct_for_template(self):



