#!/usr/bin/env python
"""
Author: Sofia Rest
"""

import json, requests, re, time, sys, os
from collections import OrderedDict
from astropy.time import Time
import numpy as np
from pdastro import pdastrostatsclass, AandB, AnotB

class atlas_lc:
	def __init__(self, tnsname=None, is_averaged=False, mjd_bin_size=None, discdate=None, ra=None, dec=None):
		#self.pdastro = pdastrostatsclass()
		self.lcs = {}

		self.tnsname = tnsname
		self.is_averaged = is_averaged
		self.mjd_bin_size = mjd_bin_size
		self.discdate = discdate
		self.ra = ra
		self.dec = dec

		self.corrected_baseline_ix = None
		self.during_sn_ix = None

		self.dflux_colnames = ['duJy']

	def __str__(self):
		res = f'SN {self.tnsname} light curve'
		if self.is_averaged:
			res += f' (averaged with MJD bin size {self.mjd_bin_size})'
		res += f': RA: {self.ra}, Dec: {self.dec}, discovery date: {self.discdate}'
		return res

	# get RA, Dec, and discovery date information from TNS
	def get_tns_data(self, api_key, tns_id, bot_name):
		print(f'Obtaining RA, Dec, and/or discovery date from TNS (TNS ID: {tns_id}; TNS bot name: "{bot_name}")...')
		if tns_id == "None" or bot_name == "None":
			raise RuntimeError("# ERROR: Cannot query TNS without TNS ID and bot name! Specify these parameters in params.ini")

		try:
			url = 'https://www.wis-tns.org/api/get/object'
			json_file = OrderedDict([("objname",self.tnsname), ("objid",""), ("photometry","1"), ("spectra","1")])
			data = {'api_key':api_key,'data':json.dumps(json_file)}
			response = requests.post(url, data=data, headers={'User-Agent': 'tns_marker{"tns_id":"%s","type": "bot", "name":"%s"}' % (tns_id, bot_name)})
			json_data = json.loads(response.text,object_pairs_hook=OrderedDict)
			
			if self.ra is None:
				self.ra = json_data['data']['reply']['ra']
			if self.dec is None: 
				self.dec = json_data['data']['reply']['dec']
			discoverydate = json_data['data']['reply']['discoverydate']
		except Exception as e:
			print(json_data['data']['reply'])
			raise RuntimeError('# ERROR in get_tns_data(): '+str(e))

		date = list(discoverydate.partition(' '))[0]
		time = list(discoverydate.partition(' '))[2]
		dateobjects = Time(date+"T"+time, format='isot', scale='utc')
		if self.discdate is None:
			self.discdate = dateobjects.mjd - 20 # make sure no SN flux before discovery date in baseline indices

	# get baseline indices (any indices before the SN discovery date)
	def get_baseline_ix(self):
		if self.discdate is None:
			raise RuntimeError('ERROR: Cannot get baseline indices because discovery date is None!')
		return self.lcs[0].ix_inrange(colnames=['MJD'],uplim=self.discdate-20,exclude_uplim=True)

	# get a light curve filename for saving/loading
	def get_filename(self, filt, control_index, directory):
		# SN light curve: 				DIRECTORY/2022xxx/2022xxx.o.lc.txt
		# averaged light curve: 		DIRECTORY/2022xxx/2022xxx.o.1.00days.lc.txt
		# control light curve: 			DIRECTORY/2022xxx/controls/2022xxx_i001.o.lc.txt
		# averaged control light curve: DIRECTORY/2022xxx/controls/2022xxx_i001.o.1.00days.lc.txt

		filename = f'{directory}/{self.tnsname}'
		if control_index != 0:
			filename += '/controls'
		filename += f'/{self.tnsname}'
		if control_index != 0:
			filename += f'_i{control_index:03d}'
		filename += f'.{filt}'
		if self.is_averaged:
			filename += f'.{self.mjd_bin_size:0.2f}days'
		filename += '.lc.txt'
		
		#print(f'# Filename: {filename}')
		return filename

	# save a single light curve
	def save_lc(self, output_dir, control_index=0, filt=None, overwrite=True, keep_empty_bins=True):
		if filt is None:
			# split lc up by filt and save to two separate files
			for filt_ in ['o','c']:
				# if not keeping empty bins in averaged lc, remove all null rows; else keep all
				if self.is_averaged and not keep_empty_bins:
					ix = AandB(self.lcs[control_index].ix_not_null(colnames=['uJy']), self.lcs[control_index].ix_equal(colnames=['F'],val=filt_))
				else: 
					ix = self.lcs[control_index].ix_equal(colnames=['F'],val=filt_)
				
				self.lcs[control_index].write(filename=self.get_filename(filt_,control_index,output_dir), indices=ix, overwrite=overwrite)
		else:
			if self.is_averaged and not keep_empty_bins:
				ix = self.lcs[control_index].ix_not_null(colnames=['uJy'])
			else: 
				ix = self.lcs[control_index].getindices()
			self.lcs[control_index].write(filename=self.get_filename(filt,control_index,output_dir), indices=ix, overwrite=overwrite)

	# save SN light curve and, if necessary, control light curves
	def save(self, output_dir, filt=None, overwrite=True, keep_empty_bins=True):
		if len(self.lcs) < 1:
			print('WARNING: No light curves to save! Skipping...')
		else:
			if self.is_averaged:
				output = f'\nSaving averaged SN light curve and {len(self.lcs)-1} averaged control light curves (keep empty bins: {keep_empty_bins})...'
			else:
				output = f'\nSaving SN light curve and {len(self.lcs)-1} control light curves...'
			print(output)

			for control_index in self.lcs:
				self.save_lc(output_dir, control_index, filt=filt, overwrite=overwrite, keep_empty_bins=keep_empty_bins)

	# load a single light curve
	def load_lc(self, output_dir, filt, is_averaged=False, control_index=0):
		if (len(self.lcs) > 0) and self.is_averaged != is_averaged:
			raise RuntimeError(f'ERROR: cannot load a light curve whose is_averaged status of {is_averaged} does not match previous status of {self.is_averaged}!')
		self.is_averaged = is_averaged

		self.lcs[control_index] = pdastrostatsclass()
		self.lcs[control_index].load_spacesep(self.get_filename(filt, control_index, output_dir), delim_whitespace=True)
		self.lcs[control_index].t['Mask'] = 0

	# load SN light curve and, if necessary, control light curves for a certain filter
	def load(self, output_dir, filt, is_averaged=False, num_controls=0):
		output = f'\nLoading averaged SN light curve and {num_controls} averaged control light curves...' if self.is_averaged else f'\nLoading SN light curve and {num_controls} control light curves...'
		print(output)

		self.load_lc(output_dir, filt, is_averaged=is_averaged)
		for control_index in range(1, num_controls+1):
			self.load_lc(output_dir, filt, is_averaged=is_averaged, control_index=control_index)

		self.dflux_colnames = ['duJy'] * (num_controls+1)

	def exists(self, output_dir, filt, is_averaged=False, control_index=0):
		filename = self.get_filename(filt, control_index, output_dir)
		return os.path.isfile(filename)

	# update given indices of 'Mask' column in the light curve (SN if control index is None) with given flag(s)
	def update_mask_col(self, flag, indices, control_index=0):
		if len(indices) > 1:
			flag_arr = np.full(self.lcs[control_index].t.loc[indices,'Mask'].shape, flag)
			self.lcs[control_index].t.loc[indices,'Mask'] = np.bitwise_or(self.lcs[control_index].t.loc[indices,'Mask'].astype(int), flag_arr)
		elif len(indices) == 1:
			self.lcs[control_index].t.loc[indices,'Mask'] = int(self.lcs[control_index].t.loc[indices,'Mask']) | flag

	# get the xth percentile SN flux using given indices
	def get_xth_percentile_flux(self, percentile, indices=None):
		if indices is None or len(indices)==0:
			indices = self.lcs[0].getindices()
		return np.percentile(self.lcs[0].t.loc[indices, 'uJy'], percentile)

	def get_filt_lens(self, control_index=0):
		o_len = len(self.lcs[control_index].ix_equal(colnames=['F'],val='o'))
		c_len = len(self.lcs[control_index].ix_equal(colnames=['F'],val='c'))
		return o_len, c_len
