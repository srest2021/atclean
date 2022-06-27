#!/usr/bin/env python
"""
Author: Sofia Rest
"""

import json, requests, re, time, sys
from collections import OrderedDict
from astropy.time import Time
from pdastro import pdastrostatsclass, AorB, AnotB

class light_curve:
	def __init__(self, tnsname=None, is_averaged=False, mjdbinsize=None, discdate=None, ra=None, dec=None):
		self.tnsname = tnsname
		self.is_averaged = is_averaged
		self.mjdbinsize = mjdbinsize
		self.discdate = discdate
		self.ra = ra
		self.dec = dec
		self.pdastro = pdastrostatsclass()
		self.lcs = {}

	# get RA, Dec, and discovery date information from TNS
	def get_tns_data(self, api_key):
		print('Obtaining RA, Dec, and discovery date from TNS...')

		try:
			url = 'https://www.wis-tns.org/api/get/object'
			json_file = OrderedDict([("objname",self.tnsname), ("objid",""), ("photometry","1"), ("spectra","1")])
			data = {'api_key':api_key,'data':json.dumps(json_file)}
			response = requests.post(url, data=data, headers={'User-Agent':'tns_marker{"tns_id":104739,"type": "bot", "name":"Name and Redshift Retriever"}'})
			json_data = json.loads(response.text,object_pairs_hook=OrderedDict)
		except Exception as e:
			raise RuntimeError('ERROR in get_tns_data(): '+str(e))

		self.ra = json_data['data']['reply']['ra']
		self.dec = json_data['data']['reply']['dec']

		discoverydate = json_data['data']['reply']['discoverydate']
		date = list(discoverydate.partition(' '))[0]
		time = list(discoverydate.partition(' '))[2]
		dateobjects = Time(date+"T"+time, format='isot', scale='utc')
		self.discdate = dateobjects.mjd

	# get a light curve filename for saving
	def get_filename(self, filt, control_index, directory):
		if not self.is_averaged:
			filename = f'{directory}/{self.tnsname}/{self.tnsname}_i{control_index:03d}.{filt}.lc.txt'
		else:
			filename = f'{directory}/{self.tnsname}/{self.tnsname}_i{control_index:03d}.{filt}.{self.mjdbinsize:0.2f}days.lc.txt'
		print(f'# Filename: {filename}')
		return filename

	# save SN light curve and, if necessary, control light curves
	def save(self, output_dir, overwrite=True):
		print('Saving SN light curve')
		o_ix = self.pdastro.ix_equal(colnames=['F'],val='o')
		self.pdastro.write(filename=self.get_filename('o',0,output_dir), indices=o_ix, overwrite=overwrite)
		self.pdastro.write(filename=self.get_filename('c',0,output_dir), indices=AnotB(self.pdastro.getindices(),o_ix), overwrite=overwrite)

		if len(self.lcs) > 0:
			print('Saving control light curves')
			for control_index in range(1,len(self.lcs)):
				for filt in ['c','o']:
					filt_ix = self.lcs[control_index].pdastro.ix_equal(colnames=['F'],val=filt)
					self.lcs[control_index].pdastro.write(filename=self.get_filename(filt,control_index,output_dir), indices=filt_ix, overwrite=overwrite)

	# load SN light curve and, if necessary, control light curves for a certain filter
	def load(self, filt, input_dir, num_controls=None):
		print('Loading SN light curve')
		self.pdastro.load_spacesep(self.get_filename(filt,0,input_dir), delim_whitespace=True)

		if not(num_controls is None):
			print(f'Loading {num_controls} control light curves')
			for control_index in range(1,num_controls+1):
				self.lcs[control_index].pdastro.load_spacesep(self.get_filename(filt,control_index,input_dir), delim_whitespace=True)

	# add downloaded control light curve to control light curve dictionary
	def add_control_lc(self, control_lc):
		self.lcs[len(self.lcs)+1] = control_lc