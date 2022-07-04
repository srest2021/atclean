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

	def average_loop(self):
		args = self.define_args().parse_args()
		self.load_settings(args)

		for obj_index in range(0,len(args.tnsnames)):
			print(f'\nCOMMENCING AVERAGING LOOP FOR SN {args.tnsnames[obj_index]}')

if __name__ == "__main__":
	average_atlas_lc = average_atlas_lc()
	average_atlas_lc.average_loop()
