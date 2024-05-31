#!/usr/bin/env python

import json, argparse
import math
import sys
import pandas as pd
import numpy as np
from typing import Dict, List, Callable, Tuple

from pdastro import pdastrostatsclass

SIM_TABLE_REQUIRED_COLUMNS = ['model_name', 'filename']

class SimTable(pdastrostatsclass):
	def __init__(self, **kwargs):
		pdastrostatsclass.__init__(self, **kwargs)

	def add_row(self, parameters: Dict, filename=None):
		pass

	# update a certain row of the table
	# data should be structured like: data = {'column1': value1, 'column2': value2}
	def update_row_at_index(self, index, data: Dict):
		self.t.loc[index, data.keys()] = np.array(list(data.values()))

class SimTables:
	def __init__(self, peak_appmags: List):
		self.d: Dict[str, SimTable] = None
		self.peak_appmags = peak_appmags

	def setup(self, ):
		pass

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

def generate_sim_detec_tables(model_name, parsed_params, filename=None, 
															mjd_column_name=None, mag_column_name=None, flux_column_name=None):
	pass

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

	filename, mjd_column_name, mag_column_name, flux_column_name = None, False, False, False
	if not args.model_name == 'Gaussian':
		try:
			filename = model_settings['filename']
			mjd_column_name = model_settings['mjd_column_name']
			mag_column_name = model_settings['mag_column_name']
			flux_column_name = model_settings['flux_column_name']

			if not mag_column_name and not flux_column_name:
				raise RuntimeError(f'ERROR: Model must have mag or flux column.')
		except Exception as e:
			raise RuntimeError(f'ERROR: {str(e)}')
			
	generate_sim_detec_tables(args.model_name, 
													 	parsed_params, 
														filename=filename, 
														mjd_column_name=mjd_column_name, 
														mag_column_name=mag_column_name, 
														flux_column_name=flux_column_name)