#!/usr/bin/env python

from abc import ABC, abstractmethod
import argparse
from typing import Tuple
import numpy as np
import pandas as pd

from generate_sim_table import load_json_config
from lightcurve import SimDetecLightCurve, SimDetecSupernova

# convert flux to magnitude 
def flux2mag(flux:float):
	return -2.5 * np.log10(flux) + 23.9

# convert magnitude to flux
def mag2flux(mag:float):
	return 10 ** ((mag - 23.9) / -2.5)

class Simulation(ABC):
	def __init__(self, model_name=None, **kwargs):
		"""
		Initialize the Simulation object.
		"""
		self.model_name = model_name
		self.peak_mjd = None

	@abstractmethod
	def get_sim_flux(self, mjds, peak_mjd=None, t0_mjd=None, **kwargs):
		"""
		Compute the simulated flux for given MJDs.

		:param mjds: List or array of MJDs.
		:param peak_mjd: MJD where the peak apparent magnitude should occur (optional).
		:param t0_mjd: MJD where the start of the model should occur (optional).

		:return: An array of flux values corresponding to the input MJDs.
		"""
		pass

	# @abstractmethod
	# def get_row_params(self):
	# 	"""
	# 	Construct a dictionary of parameters to add to a row in a SimDetecTable.
	# 	Dictionary keys should function as column names and dictionary values 
	# 	as column values.
	# 	"""
	# 	parameters = {
	# 		'name': self.name,
	# 		'peak_mjd': self.peak_mjd
	# 	}
	# 	return parameters

	def __str__(self):
		return f'Simulation with model name {self.model_name}'
	
# TODO: documentation
class SimDetecLoop(ABC):
	def __init__(self, output_dir:str, tables_dir:str, sigma_kerns:List, **kwargs):
		self.output_dir = output_dir
		self.tables_dir = tables_dir
		self.sigma_kerns = sigma_kerns

		self.sn:SimDetecSupernova = None
		self.e:EfficiencyTable = None
		self.sd:SimDetecTables = None

	@abstractmethod
	def set_peak_mags_and_fluxes(self, peak_mag_min:float=23.0, peak_mag_max:float=16.0, n_peaks:int=20, **kwargs):
		self.peak_appmags = list(np.linspace(peak_mag_min, peak_mag_max, num=n_peaks))
		self.peak_fluxes = list(map(mag2flux, self.peak_appmags))
	
	@abstractmethod
	def load_sn(self, data_dir, tnsname, num_controls, mjdbinsize=1.0, filt='o'):
		self.sn.load_all(data_dir, num_controls=num_controls)
		self.sn.remove_rolling_sums()
		self.sn.remove_simulations()

	@abstractmethod
	def modify_light_curve(self, control_index=0):
		pass

	@abstractmethod
	def get_max_fom(self, sim_lc:SimDetecLightCurve, **kwargs) -> Tuple[float, float]:
		pass

	@abstractmethod
	def add_row_to_sd(self, sigma_kern, peak_appmag, index, sim:Simulation, control_index, max_fom, max_fom_mjd):
		pass
		# sim_params = sim.get_row_params()
		# other_params = {
		# 	'control_index': control_index,
		# 	'max_fom': max_fom,
		# 	'max_fom_mjd': max_fom_mjd
		# }
		# row = other_params.update(sim_params)
		# self.sd.update_row_at_index(sigma_kern, peak_appmag, index, row)

	@abstractmethod
	def calculate_efficiencies(self):
		pass

	@abstractmethod
	def loop(self):
		pass
	
# define command line arguments
def define_args(parser=None, usage=None, conflict_handler='resolve'):
	if parser is None:
		parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
		
	parser.add_argument('-f', '--config_file', default='simulation_settings.json', type=str, help='file name of JSON file with general settings')
	parser.add_argument('-f', '--model_config_file', default='model_settings.json', type=str, help='file name of JSON file with model settings')
	parser.add_argument('-e', '--efficiencies', default=False, action='store_true', help='calculate efficiencies using FOM limits')

	return parser

if __name__ == "__main__":
	args = define_args().parse_args()
	config = load_json_config(args.config_file)
	model_config = load_json_config(args.model_config_file)

	sn_info = config["sn_info"]
	
	data_dir = config["data_dir"]
	sim_tables_dir = config["sim_tables_dir"]
	detec_tables_dir = config["detec_tables_dir"]

	sigma_kerns = config["sigma_kerns"]
	skip_control_ix = config["skip_control_ix"]