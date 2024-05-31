#!/usr/bin/env python

import itertools
import json, argparse
import math
import sys
import pandas as pd
import numpy as np
from typing import Dict, List, Callable, Tuple

from pdastro import pdastrostatsclass

SIM_TABLE_REQUIRED_COLUMNS = ['model_name', 'filename']

# define command line arguments
def define_args(parser=None, usage=None, conflict_handler='resolve'):
	if parser is None:
		parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
	parser.add_argument("model_name", type=str, help='name of model to use')
	parser.add_argument('-f', '--config_file', default='simulation_settings.json', type=str, help='file name of JSON file with general settings')
	parser.add_argument('-f', '--model_config_file', default='model_settings.json', type=str, help='file name of JSON file with model settings')
	return parser

# load the JSON config file
def load_json_config(config_file: str):
		try:
			print(f'Loading config file at {config_file}...')
			with open(config_file) as cfg:
				return json.load(cfg)
		except Exception as e:
			raise RuntimeError(f'ERROR: Could not load config file at {config_file}: {str(e)}')

def parse_range_param(minval: float, maxval: float, step: float):
	print(f'Setting to range from {minval} to {maxval} with step size {step}')
	if maxval <= minval:
		raise RuntimeError('ERROR: maxval must be greater than minval.')
	if step > abs(maxval - minval):
		raise RuntimeError('ERROR: step size cannot be greater than the difference between minval and maxval.')
	return list(np.arange(minval, maxval, step))

def parse_logrange_param(minval: float, maxval: float, base:int, n:int, to_int=False):
	print(f'Generating {n}-length {"integer" if to_int else "float"} range using log base {base}')
	if maxval <= minval:
		raise RuntimeError('ERROR: maxval must be greater than minval.')
	minlog = np.log(minval) / np.log(base)
	maxlog = np.log(maxval) / np.log(base)
	res = list(np.logspace(minlog, maxlog, num=n, base=base))
	if to_int:
		res = [round(num) for num in res]
	return res

def parse_random_param(minval: float, maxval: float, n:int, to_int=False):
	print(f'Generating {n}-length random list')
	if maxval <= minval:
		raise RuntimeError('ERROR: maxval must be greater than minval.')
	res = list(np.random.uniform(minval, maxval, n))
	if to_int:
		res = [round(num) for num in res]
	return res

def filter_valid_draws(valid_ranges: List[List[float]], draws: List):
	return [
		draw for draw in draws
		if any(valid_range[0] <= draw <= valid_range[1] for valid_range in valid_ranges)
	]

def rec_get_valid_draws(valid_ranges: List[List[float]], n: int):
	m = 2 * n
	m_draws = list(np.random.uniform(valid_ranges[0][0], valid_ranges[-1][1], m))
	valid_draws = filter_valid_draws(valid_ranges, m_draws)
	
	m_prime = len(valid_draws)
	if m_prime >= n:
		return valid_draws
	return valid_draws + rec_get_valid_draws(valid_ranges, n - m_prime)

def parse_random_inrange_param(valid_ranges: List[List[float]], n: int):
	print(f'Generating {n}-length random list within the following valid ranges: {valid_ranges}')
	return rec_get_valid_draws(valid_ranges, n)

def parse_param(param_name: str, param_info):
	print(f'\nParsing parameter {param_name}:')
	
	if param_info['type'] == 'list':
		print(f'Setting to list')
		res = param_info['list']
	
	elif param_info['type'] == 'range':
		res = parse_range_param(param_info['range']['minval'],
														param_info['range']['maxval'],
														param_info['range']['step'])
	
	elif param_info['type'] == 'logrange':
		res = parse_logrange_param(param_info['logrange']['minval'],
														 	 param_info['logrange']['maxval'],
															 param_info["logrange"]["base"],
															 param_info['logrange']['n'],
															 to_int=param_info['logrange']['to_int'])

	elif param_info['type'] == 'random':
		res = parse_random_param(param_info['random']['minval'],
														 param_info['random']['maxval'],
														 param_info['random']['n'],
														 to_int=param_info['random']['to_int'])
	
	elif param_info['type'] == 'random_inrange':
		res = parse_random_inrange_param(param_info['random_inrange']['valid_ranges'], 
																	 	 param_info['random_inrange']['n'])

	else:
		raise RuntimeError('ERROR: Type must be one of the following: list, range, logrange, random, random_inrange.')

	print(f'Result: {res}')
	return res

class SimTable(pdastrostatsclass):
	def __init__(self, peak_appmag, **kwargs):
		pdastrostatsclass.__init__(self, **kwargs)
		self.peak_appmag = peak_appmag

	def add_row(self, data: Dict):
		self.t = pd.concat([self.t, pd.DataFrame([data])], ignore_index=True)

	# update a certain row of the table
	# data should be structured like: data = {'column1': value1, 'column2': value2}
	def update_row_at_index(self, index, data: Dict):
		self.t.loc[index, data.keys()] = np.array(list(data.values()))

	def save(self, tables_dir):
		filename = f'simtable_{self.peak_appmag:0.2f}.txt'
		print(f'Saving SimTable {filename}...')
		filename = f'{tables_dir}/{filename}'
		self.write(filename=filename, overwrite=True, index=False)

class SimTables:
	def __init__(self, peak_appmags: List, model_name: str):
		self.d: Dict[str, SimTable] = {None}
		self.peak_appmags = [round(peak_appmag,2) for peak_appmag in peak_appmags]
		self.model_name = model_name

	def generate(self, parsed_params: Dict[str, List], filename=None, mjd_colname=None, mag_colname=None, flux_colname=None):
		del parsed_params['peak_appmag']
		num_rows = 0
		for param_name in parsed_params:
			num_rows += len(parsed_params[param_name])

		row = {
			'model_name': self.model_name,
			'filename': np.nan if filename is None else filename
		}
		if not mjd_colname is False:
			row['mjd_colname'] = np.nan if mjd_colname is None else mjd_colname
		if not mag_colname is False:
			row['mag_colname'] = np.nan if mag_colname is None else mag_colname
		if not flux_colname is False:
			row['flux_colname'] = np.nan if flux_colname is None else flux_colname
		
		print()
		self.d = {}
		for peak_appmag in self.peak_appmags:
			print(f'Generating {num_rows}-length SimTable for peak_appmag={peak_appmag}...')
			self.d[peak_appmag] = SimTable(peak_appmag)
			
			combinations = list(itertools.product(*(parsed_params.values())))
			combinations_dicts = [dict(zip(parsed_params.keys(), combo)) for combo in combinations]

			for combination in combinations_dicts:
				combination.update(row)
				self.d[peak_appmag].add_row(combination)

		print('Success')

	def save(self, tables_dir):
		print(f'\nSaving SimTables in directory: {tables_dir}')
		for peak_appmag in self.peak_appmags:
			self.d[peak_appmag].save(tables_dir)
		print('Success')

if __name__ == "__main__":
	args = define_args().parse_args()
	config = load_json_config(args.config_file)
	model_config = load_json_config(args.model_config_file)

	try:
		model_settings = model_config[args.model_name]
	except Exception as e:
		raise RuntimeError(f'ERROR: Could not find model {args.model_name} in model config file: {str(e)}')

	parsed_params = {}
	for param_name in model_settings['parameters']:
		parsed_params[param_name] = parse_param(param_name, model_settings['parameters'][param_name])

	if not 'peak_appmag' in parsed_params:
		raise RuntimeError('ERROR: Parameters must include peak apparent magnitude (\"peak_appmag\").')

	filename, mjd_colname, mag_colname, flux_colname = None, False, False, False
	if not args.model_name == 'Gaussian':
		try:
			filename = model_settings['filename']
			mjd_colname = model_settings['mjd_column_name']
			mag_colname = model_settings['mag_column_name']
			flux_colname = model_settings['flux_column_name']

			if mjd_colname is False:
				raise RuntimeError(f'ERROR: Model must have an MJD column. Please set the field to null or the correct column name.')
			if mag_colname is False and flux_colname is False:
				raise RuntimeError(f'ERROR: Model must have either mag or flux column. Please set one or both fields to null or the correct column name.')
		except Exception as e:
			raise RuntimeError(f'ERROR: {str(e)}')
		
	sim_tables = SimTables(parsed_params['peak_appmag'], args.model_name)
	sim_tables.generate(parsed_params,
										 	filename=filename, 
									    mjd_colname=mjd_colname, 
									 		mag_colname=mag_colname, 
									 		flux_colname=flux_colname)
	sim_tables.save(config['sim_tables_dir'])