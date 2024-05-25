#!/usr/bin/env python

"""
This script adds series of simulated Gaussian events to the control light curves of a target SN.
We do this for different rolling sum kernel sizes ("sigma_kerns"). 
For each sigma_kern, we vary the peak magnitude ("peak_mag") and the kernel size ("sigma_sims") of the simulated event.

- Open the config file simulation_settings.json and edit as needed.
	- To only simulate events within observation seasons:
		- add MJD ranges of observation seasons as lists [STARTMJD, ENDMJD] to "observation_seasons"
		- ./simulate.py -m
	- To automatically calculate efficiencies:
		- add detection limits to "sim_settings" for each sigma_kern object
		- ./simulate.py -e
- To load a config file with a different file name: ./simulate.py --cfg_filename simulation_settings_copy.json
- After script finishes running, open simulation_analysis.ipynb and load in light curves and saved tables to get a walkthrough analysis.
"""

from typing import Dict
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
		
		# interpolate gaussian
		fn = interp1d(g[0],g[1],bounds_error=False,fill_value=0)
		fn = fn(mjds)
		return fn 
	
	def __str__(self):
		return f'Gaussian with peak app mag {self.peak_appmag:0.2f} and sigma_sim {self.sigma}'
	
"""
SIMULATED ERUPTION FROM LIGHT CURVE
"""

class Eruption:
	def __init__(self, filename, sigma=2.8):
		self.peak_appmag = None
		self.sigma = sigma

		self.t = None 
		self.load(filename)

	def load(self, filename):
		print(f'Loading eruption lc at {filename}...')

		try:
			self.t = pd.read_table(filename,delim_whitespace=True,header=None)
			self.t = self.t.rename(columns={0: "MJD", 1: "m"})
		except Exception as e:
			raise RuntimeError(f'ERROR: Could not load eruption at {filename}: {str(e)}')
	
		# app mag -> flux
		self.t['uJy'] = self.t['m'].apply(lambda mag: mag2flux(mag))

	# get interpolated function of eruption light curve
	def erup2fn(self, mjds, peak_mjd, peak_appmag):
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

class SimDetecLightCurve(AveragedLightCurve):
	def __init__(self, control_index=0, filt='o', mjdbinsize=1.0, **kwargs):
		AveragedLightCurve.__init__(self, control_index, filt, mjdbinsize, **kwargs)
		#self.sigma_kern = None

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
			mjd_ix = []
			for i in range(len(self.t)):
				if in_valid_season(self.t.loc[i,'peak_mjd'], valid_seasons):
					mjd_ix.append(i)
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
	def __init__(self, sigma_kerns, peak_appmags, peak_fluxes, sigma_sims, fom_limits=None, **kwargs):
		pdastrostatsclass.__init__(self, **kwargs)