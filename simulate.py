#!/usr/bin/env python

"""
This script adds series of simulated Gaussian events to the control light curves of a target SN.
We do this for different rolling sum kernel sizes ("sigma_kerns"). 
For each sigma_kern, we vary the peak magnitude ("peak_mag") and the kernel size ("sigma_sims") of the simulated event.

- Open the config file simulation_settings.json and edit as needed.
	- To only simulate events within observation seasons:
		- add MJD ranges of observation seasons as lists [STARTMJD, ENDMJD] to "observation_seasons"
		- ./simulate.py -m
	- To automatically calculate efficiencies using the FOM limits in the config file:
		- add detection limits to "sim_settings" for each sigma_kern object
		- ./simulate.py -e
- To load a config file with a different file name: ./simulate.py -f simulation_settings_2.json
- After script finishes running, open simulation_analysis.ipynb and load in light curves and saved tables to get a walkthrough analysis.
"""

from typing import Dict, List
import pandas as pd
import numpy as np
import sys, random, re, os, argparse, json
from copy import deepcopy
from scipy.interpolate import interp1d
from astropy.modeling.functional_models import Gaussian1D
from download import make_dir_if_not_exists
from lightcurve import AveragedSupernova, AveragedLightCurve, AandB, AnotB, AorB
from pdastro import pdastrostatsclass

"""
UTILITY
"""

# convert flux to magnitude 
def flux2mag(flux):
	return -2.5 * np.log10(flux) + 23.9

# convert magnitude to flux
def mag2flux(mag):
	return 10 ** ((mag - 23.9) / -2.5)

def get_valid_ix(table, colname, mjd_ranges):
  valid_ix = []
  for i in range(len(table)):
    if in_valid_season(table.loc[i,colname], mjd_ranges):
      valid_ix.append(i)
  return valid_ix

# check if MJD is within valid MJD season
def in_valid_season(mjd, valid_seasons):
	in_season = False
	for season in valid_seasons:
		if mjd <= season[1] and mjd >= season[0]:
			in_season = True
	return in_season

# generate list of peak fluxes and app mags
def generate_peaks(peak_mag_min, peak_mag_max, n_peaks):
	peak_mags = list(np.linspace(peak_mag_min, peak_mag_max, num=n_peaks))
	peak_fluxes = list(map(mag2flux, peak_mags))

	peak_mags = [round(item, 2) for item in peak_mags]
	peak_fluxes = [round(item, 2) for item in peak_fluxes]

	return peak_mags, peak_fluxes

"""
ASYM GAUSSIAN
Adapted from A. Rest
"""

class Gaussian:
	def __init__(self, sigma, peak_appmag):
		self.peak_appmag = peak_appmag
		self.sigma = sigma
		self.g = self.new_gaussian(mag2flux(peak_appmag), sigma)

	def new_gaussian(self, peak_flux, sigma):
		x = np.arange(-100,100,.01)
		g1 = Gaussian1D(amplitude=peak_flux, stddev=sigma)(x)
		g2 = Gaussian1D(amplitude=peak_flux, stddev=sigma)(x)

		ind = np.argmin(abs(x))
		g3 = np.copy(g1)
		g3[ind:] = g2[ind:]
		gauss = np.array([x,g3])
		return gauss
	
	# get interpolated function of gaussian at peak MJD (peak_mjd) and match to time array (mjds)
	def gauss2fn(self, mjds, peak_mjd):
		g = deepcopy(self.g)
		g[0,:] += peak_mjd
		
		fn = interp1d(g[0],g[1],bounds_error=False,fill_value=0)
		fn = fn(mjds)
		return fn 
	
	def __str__(self):
		return f'Gaussian with peak app mag {self.peak_appmag:0.2f} and sigma_sim {self.sigma}'
	
"""
SIMULATED MODEL FROM LIGHT CURVE
"""

class Model:
	def __init__(self, filename, sigma=2.8):
		self.peak_appmag = None
		self.sigma = sigma

		self.t = None 
		self.load(filename)

	def load(self, filename):
		print(f'Loading model light curve at {filename}...')

		try:
			self.t = pd.read_table(filename,delim_whitespace=True,header=None)
			self.t = self.t.rename(columns={0: "MJD", 1: "m"})
		except Exception as e:
			raise RuntimeError(f'ERROR: Could not load model at {filename}: {str(e)}')
	
		# app mag -> flux
		self.t['uJy'] = self.t['m'].apply(lambda mag: mag2flux(mag))

	# get interpolated function of eruption light curve
	def model2fn(self, mjds, peak_mjd, peak_appmag):
		self.peak_appmag = peak_appmag

		peak_idx = self.t['m'].idxmin() # get peak appmag

		self.t['MJD'] -= self.t.loc[peak_idx,'MJD'] # put peak appmag at days=0
		self.t['MJD'] += peak_mjd # put peak appmag at days=peak_mjd

		# scale
		self.t['uJy'] *= mag2flux(peak_appmag)/self.t.loc[peak_idx, 'uJy']
		
		# flux -> app mag
		self.t['m'] = self.t['uJy'].apply(lambda flux: flux2mag(flux)) 
		
		# interpolate lc
		fn = interp1d(self.t['MJD'], self.t['uJy'], bounds_error=False, fill_value=0)
		fn = fn(mjds)
		return fn

"""
ADD SIMULATIONS AND APPLY ROLLING SUM TO AVERAGED LIGHT CURVE 
"""

class SimDetecSupernova(AveragedSupernova):
	def __init__(self, tnsname:str=None, mjdbinsize:float=1.0, filt:str='o'):
		AveragedSupernova.__init__(self, tnsname=tnsname, mjdbinsize=mjdbinsize, filt=filt)
		self.avg_lcs: Dict[int, SimDetecLightCurve] = {}

	def apply_rolling_sums(self, sigma_kern:float, flag=0x800000):
		for control_index in range(self.num_controls+1):
			self.avg_lcs[control_index].apply_rolling_sum(sigma_kern, flag=flag)

	def remove_rolling_sums(self):
		for control_index in range(self.num_controls+1):
			self.avg_lcs[control_index].remove_rolling_sum()

	def remove_simulations(self):
		for control_index in range(self.num_controls+1):
			self.avg_lcs[control_index].remove_simulations()

	def add_simulations(self, control_index=0):
		self.avg_lcs[control_index].remove_simulations()
		# TODO

class SimDetecLightCurve(AveragedLightCurve):
	def __init__(self, control_index=0, filt='o', mjdbinsize=1.0, **kwargs):
		AveragedLightCurve.__init__(self, control_index, filt, mjdbinsize, **kwargs)
		self.cur_sigma_kern = None

	def get_valid_seasons_ix(self, valid_seasons, column='MJDbin'):
		return get_valid_ix(self.t, column, valid_seasons)
	
	# drop rolling sum columns
	def remove_rolling_sum(self):
		self.cur_sigma_kern = None

		dropcols = []
		for col in ['__tmp_SN','SNR','SNRsum','SNRsumnorm']:
			if col in self.t.columns:
				dropcols.append(col)
		if len(dropcols) > 0:
			self.t.drop(columns=dropcols,inplace=True)
	
	def remove_simulations(self):
		dropcols = []
		for col in ['__tmp_SN','uJysim','SNRsim','simLC','SNRsimsum']:
			if col in self.t.columns:
				dropcols.append(col)
		if len(dropcols) > 0:
			self.t.drop(columns=dropcols,inplace=True)

	def apply_rolling_sum(self, sigma_kern, indices=None, flag=0x800000, verbose=False):
		if indices is None:
			indices = self.getindices()
		if len(indices) < 1:
			raise RuntimeError('ERROR: not enough measurements to apply simulated gaussian')
		good_ix = AandB(indices, self.ix_unmasked('Mask', flag))

		self.cur_sigma_kern = sigma_kern

		self.remove_rolling_sum()
		self.t.loc[indices, 'SNR'] = 0.0
		self.t.loc[good_ix,'SNR'] = self.t.loc[good_ix,'uJy']/self.t.loc[good_ix,'duJy']

		new_gaussian_sigma = round(sigma_kern/self.mjd_bin_size)
		windowsize = int(6 * new_gaussian_sigma)
		halfwindowsize = int(windowsize * 0.5) + 1
		if verbose:
			print(f'Sigma: {sigma_kern:0.2f} days; MJD bin size: {self.mjd_bin_size:0.2f} days; sigma: {new_gaussian_sigma:0.2f} bins; window size: {windowsize} bins')

		# calculate the rolling SNR sum
		l = len(self.t.loc[indices])
		dataindices = np.array(range(l) + np.full(l, halfwindowsize))
		temp = pd.Series(np.zeros(l + 2*halfwindowsize), name='SNR', dtype=np.float64)
		temp[dataindices] = self.t.loc[indices,'SNR']
		SNRsum = temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=new_gaussian_sigma)
		self.t.loc[indices,'SNRsum'] = list(SNRsum[dataindices])
		
		# normalize it
		norm_temp = pd.Series(np.zeros(l + 2*halfwindowsize), name='norm', dtype=np.float64)
		norm_temp[np.array(range(l) + np.full(l, halfwindowsize))] = np.ones(l)
		norm_temp_sum = norm_temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=new_gaussian_sigma)
		self.t.loc[indices,'SNRsumnorm'] = list(SNRsum.loc[dataindices] / norm_temp_sum.loc[dataindices] * max(norm_temp_sum.loc[dataindices]))

	def add_simulation(self, lc, good_ix, sim_flux, cur_sigma_kern=None, verbose=False):
		if cur_sigma_kern is None:
			cur_sigma_kern = self.cur_sigma_kern
		if cur_sigma_kern is None:
			raise RuntimeError('ERROR: No current sigma kern passed as argument or stored during previously applied rolling sum.')
		
		lc.t.loc[good_ix,'uJysim'] += sim_flux 

		# make sure all bad rows have SNRsim = 0.0 so they have no impact on the rolling SNRsum
		lc.t['SNRsim'] = 0.0
		# include only simulated flux in the SNR
		lc.t.loc[good_ix,'SNRsim'] = lc.t.loc[good_ix,'uJysim']/lc.t.loc[good_ix,'duJy']

		new_gaussian_sigma = round(cur_sigma_kern/self.mjdbinsize)
		windowsize = int(6 * new_gaussian_sigma)
		halfwindowsize = int(windowsize * 0.5) + 1
		if verbose:
			print(f'Sigma: {cur_sigma_kern:0.2f} days; MJD bin size: {self.mjdbinsize:0.2f} days; new sigma: {new_gaussian_sigma:0.2f} bins; window size: {windowsize} bins')

		# calculate the rolling SNR sum for SNR with simulated flux
		l = len(self.t)
		dataindices = np.array(range(l) + np.full(l, halfwindowsize))
		temp = pd.Series(np.zeros(l + 2*halfwindowsize), name='SNRsim', dtype=np.float64)
		temp[dataindices] = lc.t['SNRsim']
		SNRsimsum = temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=new_gaussian_sigma)
		lc.t['SNRsimsum'] = list(SNRsimsum.loc[dataindices])

		return lc

	def add_model(self, model:Model, peak_mjd, peak_appmag, cur_sigma_kern=None, flag=0x800000, verbose=False):
		if verbose:
			print(f'Adding simulated eruption: peak MJD = {peak_mjd:0.2f} MJD; peak app mag = {peak_appmag:0.2f}; sigma = {model.sigma:0.2f} days')

		lc = deepcopy(self)
		good_ix = AandB(lc.getindices(), lc.ix_unmasked('Mask', flag))
		sim_flux = model.model2fn(lc.t.loc[good_ix,'MJD'], peak_mjd, peak_appmag)

		return self.add_simulation(lc, good_ix, sim_flux, cur_sigma_kern=cur_sigma_kern, verbose=verbose)

	def add_gaussian(self, gaussian:Gaussian, peak_mjd, cur_sigma_kern=None, flag=0x800000, verbose=False):
		if verbose:
			print(f'Adding simulated gaussian: peak MJD = {peak_mjd:0.2f} MJD; peak app mag = {gaussian.peak_appmag:0.2f}; sigma = {gaussian.sigma:0.2f} days')

		lc = deepcopy(self)
		good_ix = AandB(lc.getindices(), lc.ix_unmasked('Mask', flag))
		sim_flux = gaussian.gauss2fn(lc.t.loc[good_ix,'MJD'], peak_mjd)
			
		return self.add_simulation(lc, good_ix, sim_flux, cur_sigma_kern=cur_sigma_kern, verbose=verbose)
		
"""
SIMULATION DETECTION AND EFFICIENCY TABLES
"""

# get simulation detection dictionary (sd) key for a sigma_kern peak_appmag pair
def get_key(sigma_kern, peak_appmag=None, peak_flux=None):
	if peak_appmag is None and peak_flux is None:
		raise RuntimeError('ERROR: Cannot get the SimDetecTable key without peak app mag or flux.')
	if not peak_flux is None:
		peak_appmag = flux2mag(peak_flux)
	return f'{sigma_kern}_{peak_appmag:0.2f}'

# simulation detection table for a sigma_kern,peak_appmag pair
class SimDetecTable(pdastrostatsclass):
	def __init__(self, sigma_kern, num_iterations=None, peak_appmag=None, peak_flux=None, **kwargs):
		pdastrostatsclass.__init__(self, **kwargs)
		self.sigma_kern = sigma_kern
		
		if peak_appmag is None and peak_flux is None:
			raise RuntimeError('ERROR: Either peak app mag or flux required to construct SimDetecTable.')
		elif peak_appmag is None:
			peak_appmag = flux2mag(peak_flux)
		elif peak_flux is None:
			peak_flux = mag2flux(peak_appmag)
		self.peak_appmag = peak_appmag
		self.peak_flux = peak_flux

		if not num_iterations is None:
			self.setup(num_iterations)

	def setup(self, num_iterations):
		self.t = pd.DataFrame(columns=['sigma_kern', 'peak_appmag', 'peak_flux', 'peak_mjd', 'sigma_sim', 'max_fom', 'max_fom_mjd'])
		self.t['sigma_kern'] = np.full(num_iterations, self.sigma_kern)
		self.t['peak_appmag'] = np.full(num_iterations, self.peak_appmag)
		self.t['peak_flux'] = np.full(num_iterations, self.peak_flux)

	def get_filename(self, tables_dir):
		key = get_key(self.sigma_kern, peak_appmag=self.peak_appmag)
		filename = f'{tables_dir}/simdetec_{key}.txt'
		return filename
	
	def load(self, tables_dir):
		filename = self.get_filename(tables_dir)
		try:
			self.load_spacesep(filename, delim_whitespace=True)
		except Exception as e:
			raise RuntimeError(f'ERROR: Could not load SimDetecTable at {filename}: {str(e)}')

	def save(self, tables_dir, overwrite=True):
		filename = self.get_filename(tables_dir)
		make_dir_if_not_exists(filename)
		self.write(filename, overwrite=overwrite)

	def get_efficiency(self, fom_limit, valid_seasons=None, sigma_sim=None):
		if sigma_sim is None: # get efficiency for all sigma_sims (entire table)
			sigma_sim_ix = self.getindices()
		else:
			sigma_sim_ix = self.ix_equal('sigma_sim', sigma_sim) #np.where(self.t['sigma_sim'] == sigma_sim)[0]
		if len(sigma_sim_ix) < 1:
			print(f'WARNING: No sim sigma matches for {sigma_sim}; returning NaN...')
			return np.nan
		
		if not valid_seasons is None:
			# get efficiency for simulations only occurring within mjd_ranges
			mjd_ix = get_valid_ix(self.t, 'peak_mjd', valid_seasons)
			sigma_sim_ix = AandB(sigma_sim_ix, mjd_ix)
		
		detected_ix = self.ix_inrange('max_fom', lowlim=fom_limit, indices=sigma_sim_ix)
		efficiency = len(detected_ix)/len(sigma_sim_ix) * 100
		return efficiency

class SimDetecTables:
	def __init__(self, sigma_kerns, num_iterations=None, peak_appmags=None, peak_fluxes=None):
		self.d:Dict[str, SimDetecTable] = None

		self.sigma_kerns = sigma_kerns

		if peak_appmags is None and peak_fluxes is None:
			raise RuntimeError('ERROR: either peak app mags or peak fluxes required to construct SimDetecTables.')
		elif peak_appmags is None:
			peak_appmags = list(map(flux2mag, peak_fluxes))
		elif peak_fluxes is None:
			peak_fluxes = list(map(mag2flux, peak_appmags))
		self.peak_appmags = peak_appmags
		self.peak_fluxes = peak_fluxes

		if not num_iterations is None:
			self.setup(num_iterations)

	def setup(self, num_iterations):
		self.d = {}
		for sigma_kern in self.sigma_kerns:
			for peak_appmag in self.peak_appmags:
				key = get_key(sigma_kern, peak_appmag=peak_appmag)
				self.d[key] = SimDetecTable(sigma_kern, 
																		peak_appmag=peak_appmag, 
																		num_iterations=num_iterations)
				
	def get_efficiency(self, sigma_kern, peak_appmag, fom_limit, valid_seasons=None, sigma_sim=None):
		key = get_key(sigma_kern, peak_appmag=peak_appmag)
		return self.d[key].get_efficiency(fom_limit, valid_seasons=valid_seasons, sigma_sim=sigma_sim)
	
	def load_all(self, tables_dir):
		self.d = {}
		for sigma_kern in self.sigma_kerns:
			for peak_appmag in self.peak_appmags:
				key = get_key(sigma_kern, peak_appmag=peak_appmag)
				self.d[key] = SimDetecTable(sigma_kern, peak_appmag=peak_appmag)
				self.d[key].load(tables_dir)

	def save_all(self, tables_dir):
		for table in self.d.values():
			table.save(tables_dir)

class EfficiencyTable(pdastrostatsclass):
	def __init__(self, sigma_kerns, sigma_sims, peak_appmags=None, peak_fluxes=None, fom_limits=None, **kwargs):
		pdastrostatsclass.__init__(self, **kwargs)

		self.sigma_kerns = sigma_kerns
		self.set_sigma_sims(sigma_sims)
		self.set_fom_limits(fom_limits)

		if peak_appmags is None and peak_fluxes is None:
			raise RuntimeError('ERROR: either peak app mags or peak fluxes required to construct an EfficiencyTable.')
		elif peak_appmags is None:
			peak_appmags = list(map(flux2mag, peak_fluxes))
		elif peak_fluxes is None:
			peak_fluxes = list(map(mag2flux, peak_appmags))
		self.peak_appmags = peak_appmags
		self.peak_fluxes = peak_fluxes

		self.setup()

	def setup(self):
		self.t = pd.DataFrame(columns=['sigma_kern', 'peak_appmag', 'peak_flux', 'sigma_sim'])

		for sigma_kern in self.sigma_kerns: 
			n = len(self.peak_appmags) * len(self.sigma_sims[sigma_kern])
			
			df = pd.DataFrame(columns=['sigma_kern', 'peak_appmag', 'peak_flux', 'sigma_sim'])
			df['sigma_kern'] = np.full(n, sigma_kern)
			df['peak_appmag'] = np.repeat(self.peak_appmags, len(self.sigma_sims[sigma_kern]))
			df['peak_flux'] = np.repeat(self.peak_fluxes, len(self.sigma_sims[sigma_kern]))
			
			j = 0
			while(j < n):
				for sigma_sim in self.sigma_sims[sigma_kern]:
					df.loc[j, 'sigma_sim'] = sigma_sim
					j += 1
			self.t = pd.concat([self.t, df], ignore_index=True)

	def set_sigma_sims(self, sigma_sims):	
		if isinstance(sigma_sims, List):
			self.sigma_sims = {}
			for i in range(len(self.sigma_kerns)):
				self.sigma_sims[self.sigma_kerns[i]] = sigma_sims[i]
		else:
			if len(sigma_sims) != len(self.sigma_kerns):
				raise RuntimeError('ERROR: Each entry in sigma_kerns must have a matching list in sigma_sims')
			self.sigma_sims = sigma_sims

	def set_fom_limits(self, fom_limits):
		if fom_limits is None:
			return
				
		if isinstance(fom_limits, list):
			self.fom_limits = {}
			for i in range(len(self.sigma_kerns)):
				self.fom_limits[self.sigma_kerns[i]] = fom_limits[i]
		else:
			if len(fom_limits) != len(self.sigma_kerns):
				raise RuntimeError('ERROR: Each entry in sigma_kerns must have a matching list in fom_limits')
			self.fom_limits = fom_limits

	def reset_table(self):
		for col in self.t.columns:
			if re.search('^pct_detec_',col):
				self.t.drop(col, axis=1, inplace=True)

		for i in range(len(self.sigma_kerns)):
			fom_limits = self.fom_limits[self.sigma_kerns[i]]
			for fom_limit in fom_limits:
				self.t[f'pct_detec_{fom_limit:0.2f}'] = np.full(len(self.t), np.nan)
	
	def get_efficiencies(self, sd:SimDetecTables, valid_seasons=None, verbose=False):
		if self.fom_limits is None:
			raise RuntimeError('ERROR: fom_limits is None')
		
		for i in range(len(self.t)):
			sigma_kern = self.t.loc[i,'sigma_kern']
			peak_appmag = self.t.loc[i,'peak_appmag']
			sigma_sim = self.t.loc[i,'sigma_sim']

			if verbose:
				print(f'Getting efficiencies for sigma_kern {sigma_kern}, sigma_sim {sigma_sim}, peak_mag {peak_appmag}...')

			for fom_limit in self.fom_limits[sigma_kern]:
				self.t.loc[i,f'pct_detec_{fom_limit:0.2f}'] = \
					sd.get_efficiency(sigma_kern, 
											 			peak_appmag,
														fom_limit, 
														valid_seasons=valid_seasons, 
														sigma_sim=sigma_sim)
				
	def get_subset(self, sigma_kern=None, fom_limit=None, sigma_sim=None):
		colnames = ['sigma_kern', 'peak_appmag', 'peak_flux', 'sigma_sim']
		
		if sigma_kern is None:
			ix = self.getindices()
		else:
			ix = self.ix_equal('sigma_kern', sigma_kern) 

		if not sigma_sim is None: 
			ix = self.ix_equal('sigma_sim', sigma_sim, indices=ix)

		if not fom_limit is None:
			try: 
				colnames.append(f'pct_detec_{fom_limit:0.2f}')
			except Exception as e:
				raise RuntimeError(f'ERROR: no matching FOM limit {fom_limit}: {str(e)}')
		else:
			for col in self.t.columns:
				if re.search('^pct_detec_',col):
					colnames.append(col)
		
		return self.t.loc[ix,colnames]
	
	def merge_tables(self, other):
		if not isinstance(other, EfficiencyTable):
			raise RuntimeError(f'ERROR: Cannot merge EfficiencyTable with object type: {type(other)}')
		
		self.sigma_kerns += other.sigma_kerns
		self.sigma_sims += other.sigma_sims
		if not self.fom_limits is None:
			self.fom_limits.update(other.fom_limits)

		self.t = pd.concat([self.t, other.t], ignore_index=True)

	def load(self, tables_dir, filename='efficiencies.txt'):
		print(f'Loading efficiency table {filename}...')
		filename = f'{tables_dir}/{filename}'
		try:
			#self.t = pd.read_table(filename, delim_whitespace=True)
			self.load_spacesep(filename, delim_whitespace=True)
		except Exception as e:
			raise RuntimeError(f'ERROR: Could not load efficiency table at {filename}: {str(e)}')
	
	def save(self, tables_dir, filename='efficiencies.txt'):
		print(f'Saving efficiency table {filename}...')
		filename = f'{tables_dir}/{filename}'
		#self.t.to_string(filename, index=False)
		self.write(filename=filename, overwrite=True, index=False)
	
	def __str__(self):
		return self.t.to_string()
	
"""
GENERATE AND SAVE SIMULATED DETECTION AND EFFICIENCY TABLES
"""

# define command line arguments
def define_args(parser=None, usage=None, conflict_handler='resolve'):
	if parser is None:
		parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
		
	parser.add_argument('-f', '--config_file', default='simulation_settings.json', type=str, help='file name of JSON file with settings for this class')
	parser.add_argument('-e', '--efficiencies', default=False, action='store_true', help='calculate efficiencies using FOM limits')
	parser.add_argument('-m', '--obs_seasons', default=False, action='store_true', help='only simulate events with peak MJDs within the observation season MJD ranges')

	return parser

def load_json_config(config_file):
		with open(config_file) as cfg:
			return json.load(cfg)
		
class SimDetecLoop:
	def __init__(self, output_dir:str, tables_dir:str, sigma_kerns:List, sigma_sims, 
							 peak_mag_min:float=23.0, peak_mag_max:float=16.0, n_peaks:int=20, fom_limits:Dict=None, num_iterations:int=50000, calc_efficiencies:bool=False):
		self.output_dir = output_dir
		self.tables_dir = tables_dir
		self.num_iterations = num_iterations
		self.calc_efficiencies = calc_efficiencies

		self.sn:SimDetecSupernova = None
		
		self.sigma_kerns = sigma_kerns
		self.sigma_sims = sigma_sims
		self.fom_limits = fom_limits

		self.generate_peaks(peak_mag_min, peak_mag_max, n_peaks)

		self.e:EfficiencyTable = EfficiencyTable(sigma_kerns, sigma_sims, peak_appmags=self.peak_appmags, peak_fluxes=self.peak_fluxes, fom_limits=fom_limits)
		self.sd:SimDetecTables = SimDetecTables(sigma_kerns, num_iterations=num_iterations, peak_appmags=self.peak_appmags, peak_fluxes=self.peak_fluxes)

	def generate_peaks(self, peak_mag_min, peak_mag_max, n_peaks):
		peak_mags = list(np.linspace(peak_mag_min, peak_mag_max, num=n_peaks))
		peak_fluxes = list(map(mag2flux, peak_mags))
		self.peak_appmags = [round(item, 2) for item in peak_mags]
		self.peak_fluxes = [round(item, 2) for item in peak_fluxes]

	def loop(self, tnsname, num_controls, valid_control_ix, filt='o', mjdbinsize=1.0, valid_seasons=None):
		# load SN and control light curves
		self.sn = SimDetecSupernova(tnsname=tnsname, mjdbinsize=mjdbinsize, filt=filt)
		self.sn.load_all(output_dir, num_controls=num_controls)

		for sigma_kern in self.sigma_kerns:
			print(f'\nUsing rolling sum kernel size sigma_kern={sigma_kern} days...')

			self.sn.apply_rolling_sums(sigma_kern)

if __name__ == "__main__":
	args = define_args().parse_args()
	config = load_json_config(args.config_file)

	output_dir = config["data_dir"]
	tables_dir = config["tables_dir"]
	make_dir_if_not_exists(tables_dir)
	print(f'Data directory containing SN and control light curves: {output_dir}')
	print(f'Directory to store tables in: {tables_dir}')
	print(f'\nEfficiency calculation: {args.efficiencies}')
	num_iterations = int(config["num_iterations"])
	print(f'\nNumber of simulated events per peak magnitude: {num_iterations}')

	if args.obs_seasons:
		if len(config["observation_seasons"]) < 1:
			raise RuntimeError("ERROR: Please fill out observation seasons in the config file before using the -m argument.")
		valid_seasons = config["observation_seasons"]
		print(f'\nSimulating events only within the following observation seasons: \n{valid_seasons}')
	else:
		valid_seasons = None
	
	num_controls = int(config["sn_settings"]["num_controls"])
	valid_control_ix = [i for i in range(1, num_controls+1) if not i in config["skip_control_ix"]]
	print(f'\nSimulating events only within the following control light curves: \n', valid_control_ix)

	sigma_kerns = []
	sigma_sims = {}
	fom_limits = None
	if args.efficiencies:
		fom_limits = {}
	for obj in config['sim_settings']:
		sigma_kerns.append(obj['sigma_kern'])
		sigma_sims[obj['sigma_kern']] = obj['sigma_sims']
		if args.efficiencies:
			if  len(obj['fom_limits']) < 1:
				raise RuntimeError(f'ERROR: Efficiency calculation set to {args.efficiencies}, but no FOM limits provided for sigma_kern={obj["sigma_kern"]}.')
			fom_limits[obj['sigma_kern']] = obj['fom_limits']
	print(f'\nRolling sum kernel sizes (sigma_kerns): \n{sigma_kerns}')
	print(f'\nSimulated event kernel sizes (sigma_sims) for each sigma_kern: \n{sigma_sims}')
	if args.efficiencies:
		print(f'\nFOM limits for each sigma_kern: \n{fom_limits}')

	simdetec = SimDetecLoop(output_dir, 
												 	tables_dir,
													sigma_kerns, 
													sigma_sims,
													peak_mag_min=config["peak_mag_settings"]["peak_mag_min"],
													peak_mag_max=config["peak_mag_settings"]["peak_mag_max"],
													n_peaks=int(config["peak_mag_settings"]["n_peaks"]),
													fom_limits=fom_limits,
													num_iterations=num_iterations, 
													calc_efficiencies=args.efficiencies)
	
	print(f'\nSimulation peak magnitudes: \n{simdetec.peak_appmags}')
	print(f'Simulation peak fluxes (uJy): \n{simdetec.peak_fluxes}')

	simdetec.loop(config["sn_settings"]["tnsname"],
							 	num_controls,
								valid_control_ix,
								filt=config["sn_settings"]["filt"],
								mjdbinsize=config["sn_settings"]["mjd_bin_size"],
								valid_seasons=valid_seasons)
	
