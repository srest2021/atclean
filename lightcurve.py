#!/usr/bin/env python

from typing import Dict, Type, Any
import re, json, requests, time, sys, io, copy
from astropy import units as u
from astropy.coordinates import Angle
from astropy.time import Time
from collections import OrderedDict
from pdastro import pdastrostatsclass
import numpy as np
import pandas as pd

# number of days to subtract from TNS discovery date to make sure no SN flux before discovery date
DISC_DATE_BUFFER = 20 

# required light curve column names for the script to work
REQUIRED_COLUMN_NAMES = ['MJD', 'uJy', 'duJy']

ATLAS_FILTERS = ['c', 'o']

"""
UTILITY
"""

def AandB(A,B):
	return np.intersect1d(A,B,assume_unique=False)

def AnotB(A,B):
	return np.setdiff1d(A,B)

def AorB(A,B):
	return np.union1d(A,B)

def not_AandB(A,B):
	return np.setxor1d(A,B)

class RA:
	def __init__(self, string=None):
		self.angle = None
		if string:
			self.set_angle(string)

	def set_angle(self, string):
		s = re.compile('\:')
		if isinstance(string, str) and s.search(string):
			A = Angle(string, u.hour)
		else:
			A = Angle(string, u.degree)
		self.angle:Angle = A
	
class Dec:
	def __init__(self, string=None):
		self.angle = None
		if string:
			self.set_angle(string)

	def set_angle(self, string):
		self.angle:Angle = Angle(string, u.degree)
	
class Coordinates:
	def __init__(self, ra:str|None=None, dec:str|None=None):
		self.ra:RA = RA(ra)
		self.dec:Dec = Dec(dec)

	def set_RA(self, ra):
		self.ra = RA(ra)
	
	def set_Dec(self, dec):
		self.dec = Dec(dec)

	def is_empty(self) -> bool:
		return self.ra.angle is None or self.dec.angle is None

	def __str__(self):
		if self.is_empty():
			raise RuntimeError(f'ERROR: Coordinates are empty and cannot be printed.')
		return f'RA {self.ra.angle.degree:0.14f}, Dec {self.dec.angle.degree:0.14f}'
	
def get_filename(output_dir, tnsname, filt='o', control_index=0, mjdbinsize=None, cleaned=False):
	filename = f'{output_dir}/{tnsname}'
	
	if control_index != 0:
		filename += '/controls'
	
	filename += f'/{tnsname}'
	
	if control_index != 0:
		filename += f'_i{control_index:03d}'
	
	filename += f'.{filt}'
	
	if mjdbinsize:
		filename += f'.{mjdbinsize:0.2f}days'

	if cleaned:
		filename += f'.clean'
		
	filename += '.lc.txt'
	return filename

def query_tns(tnsname, api_key, tns_id, bot_name):
	if tns_id == 'None' or bot_name == 'None':
			raise RuntimeError('ERROR: Cannot query TNS without TNS ID and bot name. Please specify these parameters in settings.ini.')
	
	try:
		url = 'https://www.wis-tns.org/api/get/object'
		json_file = OrderedDict([("objname",tnsname), ("objid",""), ("photometry","1"), ("spectra","1")])
		data = { 'api_key':api_key, 'data':json.dumps(json_file) }
		response = requests.post(url, data=data, headers={'User-Agent': 'tns_marker{"tns_id":"%s","type": "bot", "name":"%s"}' % (tns_id, bot_name)})
		json_data = json.loads(response.text,object_pairs_hook=OrderedDict)
		return json_data
	except Exception as e:
		print(json_data['data']['reply'])
		raise RuntimeError('ERROR in query_tns(): '+str(e))

def query_atlas(headers, ra, dec, min_mjd, max_mjd):
	baseurl = 'https://fallingstar-data.com/forcedphot'
	task_url = None
	while not task_url:
		with requests.Session() as s:
			resp = s.post(f"{baseurl}/queue/",headers=headers,data={'ra':ra,'dec':dec,'send_email':False,"mjd_min":min_mjd,"mjd_max":max_mjd})
			if resp.status_code == 201: 
				task_url = resp.json()['url']
				print(f'Task url: {task_url}')
			elif resp.status_code == 429:
				message = resp.json()["detail"]
				print(f'{resp.status_code} {message}')
				t_sec = re.findall(r'available in (\d+) seconds', message)
				t_min = re.findall(r'available in (\d+) minutes', message)
				if t_sec:
					waittime = int(t_sec[0])
				elif t_min:
					waittime = int(t_min[0]) * 60
				else:
					waittime = 10
				print(f'Waiting {waittime} seconds')
				time.sleep(waittime)
			else:
				print(f'ERROR {resp.status_code}')
				print(resp.text)
				sys.exit()
	
	result_url = None
	taskstarted_printed = False
	
	print('Waiting for job to start...')
	while not result_url:
		with requests.Session() as s:
			resp = s.get(task_url, headers=headers)
			if resp.status_code == 200: 
				if not(resp.json()['finishtimestamp'] is None):
					result_url = resp.json()['result_url']
					print(f"Task is complete with results available at {result_url}")
					break
				elif resp.json()['starttimestamp']:
					if not taskstarted_printed:
						print(f"Task is running (started at {resp.json()['starttimestamp']})")
						taskstarted_printed = True
					time.sleep(2)
				else:
					#print(f"Waiting for job to start (queued at {resp.json()['timestamp']})")
					time.sleep(4)
			else:
				print(f'ERROR {resp.status_code}')
				print(resp.text)
				sys.exit()
		
	with requests.Session() as s:
		if result_url is None:
			print('WARNING: Empty light curve (no data within this MJD range).')
			dfresult = pd.DataFrame(columns=['MJD','m','dm','uJy','duJy','F','err','chi/N','RA','Dec','x','y','maj','min','phi','apfit','Sky','ZP','Obs','Mask'])
		else:
			result = s.get(result_url, headers=headers).text
			dfresult = pd.read_csv(io.StringIO(result.replace("###", "")), delim_whitespace=True)
	
	return dfresult

# input/output table containing TNS names, RA, Dec, and MJD0
# (TODO: if MJD0=None, consider entire light curve as pre-SN light curve)
class SnInfoTable():
	def __init__(self, directory, filename=None):
		if filename is None:
			self.filename = f'{directory}/sninfo.txt'
		else:
			self.filename = f'{directory}/{filename}'

		try:
			print(f'Loading SN info table at {self.filename}...')
			self.t = pd.read_table(self.filename, delim_whitespace=True)
			if not 'tnsname' in self.t.columns:
				raise RuntimeError('ERROR: SN info table must have a "tnsname" column.')
			print('Success')
		except Exception as e:
			print(f'No existing SN info table at that path; creating blank table...')
			self.t = pd.DataFrame(columns=['tnsname', 'ra', 'dec', 'mjd0']) #, 'closebright_ra', 'closebright_dec'])

	def get_index(self, tnsname):
		if self.t.empty:
			return -1
		
		matching_ix = np.where(self.t['tnsname'].eq(tnsname))[0]
		if len(matching_ix) >= 2:
			print(f'WARNING: SN info table has {len(matching_ix)} matching rows for TNS name {tnsname}. Dropping duplicate rows...')
			self.t.drop(matching_ix[1:], inplace=True)
			return matching_ix[0]
		elif len(matching_ix) == 1:
			return matching_ix[0]
		else:
			return -1

	def get_row(self, tnsname):
		if self.t.empty:
			return -1, None
		
		matching_ix = np.where(self.t['tnsname'].eq(tnsname))[0]
		if len(matching_ix) >= 2:
			print(f'WARNING: SN info table has {len(matching_ix)} matching rows for TNS name {tnsname}. Dropping duplicate rows...')
			self.t.drop(matching_ix[1:], inplace=True)
			return matching_ix[0], self.t.loc[matching_ix[0],:]
		elif len(matching_ix) == 1:
			return matching_ix[0], self.t.loc[matching_ix[0],:]
		else:
			return -1, None
		
	def add_row_info(self, tnsname, coords:Coordinates=None, mjd0=None):
		if mjd0 is None:
			mjd0 = np.nan
		#row = {'tnsname':tnsname, 'ra':coords.ra.string, 'dec':coords.dec.string, 'mjd0':mjd0}
		row = {'tnsname':tnsname, 'ra':f'{coords.ra.angle.degree:0.14f}', 'dec':f'{coords.dec.angle.degree:0.14f}', 'mjd0':mjd0}
		self.add_row(row)

	def add_row(self, row):
		if len(self.t) > 0:
			matching_ix = np.where(self.t['tnsname'].eq(row['tnsname']))[0]
			if len(matching_ix) > 1:
				raise RuntimeError(f'ERROR: SN info table has {len(matching_ix)} matching rows for TNS name {row["tnsname"]}.')
		
			if len(matching_ix) > 0:
				# update existing row
				idx = matching_ix[0]
				self.t.loc[idx,:] = row
			else:
				# new row
				self.t = pd.concat([self.t, pd.DataFrame([row])], ignore_index=True)
		else:
			# new row
			self.t = pd.concat([self.t, pd.DataFrame([row])], ignore_index=True)

	def save(self):
		print(f'Saving SN info table at {self.filename}...')
		self.t.to_string(self.filename, index=False)

	def __str__(self):
		return self.t.to_string()
	
class Cut:
	def __init__(self, 
							 column:str=None, 
							 min_value:float=None, 
							 max_value:float=None, 
							 flag:int=None,
							 params:Dict[str, Any]=None):
		self.column = column
		self.min_value = min_value
		self.max_value = max_value
		self.flag = flag
		self.params = params

	def can_apply_directly(self):
		if not self.flag or not self.column or (not self.min_value and not self.max_value):
			return False
		return True
	
	def __str__(self):
		output = ''
		if self.column:
			output += f'column={self.column} '
		if self.flag:
			output += f'flag={hex(self.flag)} '
		if self.min_value:
			output += f'min_value={self.min_value} '
		if self.max_value:
			output += f'max_value={self.max_value}'
		return output
	
"""
LIGHT CURVES
"""

class Supernova:
	def __init__(self, tnsname:str=None, ra:str=None, dec:str=None, mjd0:float=None, filt='c'):
		self.tnsname = tnsname
		self.num_controls = 0
		self.coords:Coordinates = Coordinates(ra,dec)
		self.mjd0 = mjd0 
		self.filt = filt

		self.lcs: Dict[int, Type[LightCurve]] = {}
	
	def get(self, control_index=0):
		try:
			return self.lcs[control_index].t
		except: 
			raise RuntimeError(f'ERROR: Cannot get control light curve {control_index}. Num controls set to {self.num_controls} and {len(self.lcs)} lcs in dictionary.')

	def get_tns_data(self, api_key, tns_id, bot_name):
		if self.coords.is_empty() or self.mjd0 is None:
			json_data = query_tns(self.tnsname, api_key, tns_id, bot_name)

			if self.coords.is_empty():
				self.coords = Coordinates(json_data['data']['reply']['ra'], json_data['data']['reply']['dec'])
			
			if self.mjd0 is None:
				disc_date = json_data['data']['reply']['discoverydate']
				date = list(disc_date.partition(' '))[0]
				time = list(disc_date.partition(' '))[2]
				date_object = Time(date+"T"+time, format='isot', scale='utc')
				self.mjd0 = date_object.mjd - DISC_DATE_BUFFER

	def verify_mjds(self, verbose=False):
		# sort SN lc by MJD
		self.lcs[0].t.sort_values(by=['MJD'],ignore_index=True,inplace=True)

		if self.num_controls == 0:
			return
		
		if verbose:
			print('\nMaking sure SN and control light curve MJDs match up exactly...')

		sn_sorted_mjd = self.lcs[0].t['MJD'].to_numpy()

		for control_index in range(1, self.num_controls+1):
			# sort by MJD
			self.lcs[control_index].t.sort_values(by=['MJD'],ignore_index=True,inplace=True)
			control_sorted_mjd = self.lcs[control_index].t['MJD'].to_numpy()
			
			if (len(sn_sorted_mjd) != len(control_sorted_mjd)) or not np.array_equal(sn_sorted_mjd, control_sorted_mjd):
				if verbose:
					print(f'MJDs out of agreement for control light curve {control_index}, fixing...')

				only_sn_mjd = AnotB(sn_sorted_mjd, control_sorted_mjd)
				only_control_mjd = AnotB(control_sorted_mjd, sn_sorted_mjd)

				# for the MJDs only in SN, add row with that MJD to control light curve, 
				# with all values of other columns NaN
				if len(only_sn_mjd) > 0:
					for mjd in only_sn_mjd:
						self.lcs[control_index].newrow({ 'MJD':mjd, 'Mask':0 })
				
				# remove indices of rows in control light curve for which there is no MJD in the SN lc
				if len(only_control_mjd) > 0:
					ix_to_skip = []
					for mjd in only_control_mjd:
						matching_ix = self.lcs[control_index].ix_equal('MJD',mjd)
						if len(matching_ix) != 1:
							raise RuntimeError(f'ERROR: Couldn\'t find MJD={mjd} in column MJD, but should be there!')
						ix_to_skip.extend(matching_ix)
					ix = AnotB(self.lcs[control_index].getindices(),ix_to_skip)
				else:
					ix = self.lcs[control_index].getindices()
				
				# sort again
				sorted_ix = self.lcs[control_index].ix_sort_by_cols('MJD',indices=ix)
				self.lcs[control_index].t = self.lcs[control_index].t.loc[sorted_ix]
		
			self.lcs[control_index].t.reset_index(drop=True, inplace=True)

	def prep_for_cleaning(self, verbose=False):
		for control_index in range(self.num_controls+1):
			# add blank 'Mask' column
			if verbose:
				print('Adding blank \"Mask\" column...')
			self.lcs[control_index].t['Mask'] = 0

			# remove rows with duJy=0 or uJy=NaN
			self.lcs[control_index].remove_invalid_rows(verbose=verbose)

			# calculate flux/dflux column
			self.lcs[control_index].calculate_fdf_column(verbose=verbose)

		# make sure SN and control lc MJDs match up exactly
		self.verify_mjds(verbose=verbose)

	def apply_cut(self, cut:Cut):
		if not cut.can_apply_directly():
			raise RuntimeError(f'ERROR: Cannot directly apply the following cut: {cut}')
		
		sn_percent_cut = None
		for control_index in range(self.num_controls+1):
			percent_cut = self.lcs[control_index].apply_cut(cut.column, cut.flag, min_value=cut.min_value, max_value=cut.max_value)
			if control_index == 0:
				sn_percent_cut = percent_cut
		
		return sn_percent_cut
	
	def get_uncert_est_stats(self, cut:Cut):
		def get_sigma_extra(median_dflux, stdev):
			return max(0, np.sqrt(stdev**2 - median_dflux**2))
	
		stats = pd.DataFrame(columns=['control_index', 'median_dflux', 'stdev', 'sigma_extra'])
		stats['control_index'] = list(range(1, self.num_controls+1))
		stats.set_index('control_index', inplace=True)

		for control_index in range(1, self.num_controls+1):
			dflux_clean_ix = self.lcs[control_index].ix_unmasked('Mask', maskval=cut.params['uncert_cut_flag'])
			x2_clean_ix = self.lcs[control_index].ix_inrange(colnames=['chi/N'], uplim=cut.params['temp_x2_max_value'], exclude_uplim=True)
			clean_ix = AandB(dflux_clean_ix, x2_clean_ix)

			median_dflux = self.lcs[control_index].get_median_dflux(indices=clean_ix)
			
			stdev_flux = self.lcs[control_index].get_stdev_flux(indices=clean_ix)
			if stdev_flux is None:
				print(f'WARNING: Could not get flux std dev using clean indices; retrying without preliminary chi-square cut of {cut.params["temp_x2_max_value"]}...')
				stdev_flux = self.lcs[control_index].get_stdev_flux(indices=dflux_clean_ix)
				if stdev_flux is None:
					print('WARNING: Could not get flux std dev using clean indices; retrying with all indices...')
					stdev_flux = self.lcs[control_index].get_stdev_flux(control_index=control_index)
			
			sigma_extra = get_sigma_extra(median_dflux, stdev_flux)

			stats.loc[control_index, 'median_dflux'] = median_dflux
			stats.loc[control_index, 'stdev'] = stdev_flux
			stats.loc[control_index, 'sigma_extra'] = sigma_extra

		return stats
	
	def add_noise_to_dflux(self, sigma_extra):
		for control_index in range(self.num_controls+1):
			self.lcs[control_index].add_noise_to_dflux(sigma_extra)

	def calculate_control_stats(self, uncert_flag, x2_flag):
		print('Calculating control light curve statistics...')

		len_mjd = len(self.lcs[0].t['MJD'])

		# construct arrays for control lc data
		uJy = np.full((self.num_controls, len_mjd), np.nan)
		duJy = np.full((self.num_controls, len_mjd), np.nan)
		Mask = np.full((self.num_controls, len_mjd), 0, dtype=np.int32)
		
		for control_index in range(1, self.num_controls+1):
			if len(self.lcs[control_index].t) != len_mjd or not np.array_equal(self.lcs[0].t['MJD'], self.lcs[control_index].t['MJD']):
				raise RuntimeError(f'ERROR: SN lc not equal to control lc for control_index {control_index}! Rerun or debug verify_mjds().')
			else:
				uJy[control_index-1,:] = self.lcs[control_index].t['uJy']
				duJy[control_index-1,:] = self.lcs[control_index].t[self.lcs[control_index].dflux_colname]
				Mask[control_index-1,:] = self.lcs[control_index].t['Mask']

		c2_param2columnmapping = self.lcs[0].intializecols4statparams(prefix='c2_',format4outvals='{:.2f}',skipparams=['converged','i'])

		for index in range(uJy.shape[-1]):
			pda4MJD = pdastrostatsclass()
			pda4MJD.t['uJy'] = uJy[0:,index]
			pda4MJD.t[self.lcs[0].dflux_colname] = duJy[0:,index]
			pda4MJD.t['Mask'] = np.bitwise_and(Mask[0:,index], (x2_flag|uncert_flag))
			
			pda4MJD.calcaverage_sigmacutloop('uJy',
											noisecol=self.lcs[0].dflux_colname,
											maskcol='Mask',
											maskval=(x2_flag|uncert_flag),
											verbose=1, Nsigma=3.0, median_firstiteration=True)
			self.lcs[0].statresults2table(pda4MJD.statparams, c2_param2columnmapping, destindex=index)

	def apply_controls_cut(self, cut:Cut, uncert_flag, x2_flag):
		self.calculate_control_stats(uncert_flag, x2_flag)
		self.lcs[0].t['c2_abs_stn'] = self.lcs[0].t['c2_mean'] / self.lcs[0].t['c2_mean_err']

		# flag SN measurements
		self.lcs[0].flag_by_control_stats()

		# copy over SN's control cut flags to control light curve 'Mask' columns
		flags_arr = np.full(self.lcs[0].t['Mask'].shape, 
											  (cut.flag|cut.params['questionable_flag']|cut.params['x2_flag']|cut.params['stn_flag']|cut.params['Nclip_flag']|cut.params['Ngood_flag']))
		flags_to_copy = np.bitwise_and(self.lcs[0].t['Mask'], flags_arr)
		for control_index in range(1,self.num_controls+1):
			self.lcs[control_index].copy_flags(flags_to_copy)
			
		self.drop_extra_columns()

	def drop_extra_columns(self):
		for control_index in range(self.num_controls+1):
			self.lcs[control_index].drop_extra_columns()

	def load(self, input_dir, control_index=0):
		self.lcs[control_index] = LightCurve(control_index=control_index, filt=self.filt)
		self.lcs[control_index].load(input_dir, self.tnsname)
	
	def load_all(self, input_dir, num_controls=0):
		self.num_controls = num_controls
		self.load(input_dir)
		if num_controls > 0:
			for control_index in range(1, num_controls+1):
				self.load(input_dir, control_index=control_index)

	def __str__(self):
		return f'SN {self.tnsname} at {self.coords}: MJD0 = {self.mjd0}, {self.num_controls} control light curves'

class AveragedSupernova(Supernova):
	def __init__(self, tnsname:str=None, ra:str=None, dec:str=None, mjd0:float|None = None, mjdbinsize:float=1.0):
		Supernova.__init__(self, tnsname, ra, dec, mjd0)
		self.mjdbinsize = mjdbinsize

		self.avg_lcs: Dict[int, Type[AveragedLightCurve]] = {}

	def get_avg(self, control_index:int=0):
		try:
			return self.avg_lcs[control_index].t
		except:
			raise RuntimeError(f'Cannot get averaged control light curve {control_index}. Num controls set to {self.num_controls} and {len(self.avg_lcs)} lcs in dictionary.')

	def __str__(self):
		return f'Averaged SN {self.tnsname} at {self.coords}: MJD0 = {self.mjd0}, {self.num_controls} control light curves'

# will contain measurements from both filters (o-band and c-band)
class FullLightCurve:
	def __init__(self, control_index=0, ra:str=None, dec:str=None, mjd0:float=None):
		self.t = None
		self.mjd0 = mjd0
		self.coords = Coordinates(ra,dec)
		self.control_index = control_index

	def get_tns_data(self, tnsname, api_key, tns_id, bot_name):
		if self.coords.is_empty() or self.mjd0 is None or np.isnan(self.mjd0):
			print('Querying TNS for RA, Dec, and discovery date...')
			json_data = query_tns(tnsname, api_key, tns_id, bot_name)

			if self.coords.is_empty():
				self.coords = Coordinates(json_data['data']['reply']['ra'], json_data['data']['reply']['dec'])
				print(f'Setting coordinates to TNS coordinates: {self.coords}')
			
			if self.mjd0 is None or np.isnan(self.mjd0):
				disc_date = json_data['data']['reply']['discoverydate']
				date = list(disc_date.partition(' '))[0]
				time = list(disc_date.partition(' '))[2]
				date_object = Time(date+"T"+time, format='isot', scale='utc')
				self.mjd0 = date_object.mjd - DISC_DATE_BUFFER
				print(f'Setting MJD0 to TNS discovery date minus {DISC_DATE_BUFFER}: {self.mjd0}')

	# download the full light curve from ATLAS
	def download(self, headers, lookbacktime=None, max_mjd=None):
		if lookbacktime:
			min_mjd = float(Time.now().mjd - lookbacktime)
		else:
			min_mjd = 50000.0
		
		if not max_mjd:
			max_mjd = float(Time.now().mjd)

		print(f'Downloading ATLAS light curve at {self.coords} from {min_mjd} MJD to {max_mjd} MJD...')

		if min_mjd > max_mjd:
			raise RuntimeError(f'ERROR: max MJD {max_mjd} cannot be than min MJD {min_mjd}.')
		
		while(True):
			try:
				result = query_atlas(headers, self.coords.ra.angle.degree, self.coords.dec.angle.degree, min_mjd, max_mjd)
				break
			except Exception as e:
				print('Exception caught: '+str(e))
				print('Trying again in 20 seconds! Waiting...')
				time.sleep(20)
				continue
		self.t = result

	def get_filt_lens(self):
		total_len = len(self.t)
		o_len = len(np.where(self.t['F'] == 'o')[0])
		c_len = len(np.where(self.t['F'] == 'c')[0])
		return total_len, o_len, c_len

	# divide the light curve by filter and save into separate files
	def save(self, input_dir, tnsname, overwrite=False):
		if self.t is None:
			raise RuntimeError('ERROR: Cannot save light curve that hasn\'t been downloaded yet.')

		lc = LightCurve(control_index=self.control_index)
		lc.set_df(self.t)

		# sort data by mjd
		lc.t = lc.t.sort_values(by=['MJD'],ignore_index=True)

		# remove rows with duJy=0 or uJy=NaN
		dflux_zero_ix = lc.ix_equal(colnames=['duJy'], val=0)
		flux_nan_ix = lc.ix_is_null(colnames=['uJy'])
		if len(AorB(dflux_zero_ix,flux_nan_ix)) > 0:
			print(f'Deleting {len(dflux_zero_ix) + len(flux_nan_ix)} rows with duJy=0 or uJy=NaN...')
			lc.t = lc.t.drop(AorB(dflux_zero_ix,flux_nan_ix))

		for filt in ATLAS_FILTERS:
			filename = get_filename(input_dir, tnsname, filt=filt, control_index=self.control_index)
			indices = lc.ix_equal(colnames=['F'], val=filt)
			print(f'Saving downloaded light curve with filter {filt} (length {len(indices)}) at {filename}...')
			lc.save_by_filename(filename, indices=indices, overwrite=overwrite)

	def __str__(self):
		return f'Full light curve at {self.coords}: control ID = {self.control_index}, MJD0 = {self.mjd0}'

# contains either o-band or c-band measurements only
class LightCurve(pdastrostatsclass):
	def __init__(self, control_index=0, filt='o'):
		pdastrostatsclass.__init__(self)
		self.control_index = control_index
		self.filt = filt
		self.dflux_colname = 'duJy'

	def set_df(self, t:pd.DataFrame):
		self.t = copy.deepcopy(t)

	def remove_invalid_rows(self, verbose=False):
		dflux_zero_ix = self.ix_equal(colnames=['duJy'], val=0)
		flux_nan_ix = self.ix_is_null(colnames=['uJy'])
		if len(AorB(dflux_zero_ix,flux_nan_ix)) > 0:
			if verbose:
				print(f'Deleting {len(dflux_zero_ix) + len(flux_nan_ix)} rows with duJy=0 or uJy=NaN...')
			self.t.drop(AorB(dflux_zero_ix,flux_nan_ix), inplace=True)

	def calculate_fdf_column(self, verbose=False):
		# replace infs with NaNs
		if verbose:
			print('Replacing infs with NaNs...')
		self.t.replace([np.inf, -np.inf], np.nan, inplace=True)

		# calculate flux/dflux
		if verbose:
			print('Calculating flux/dflux...')
		self.t[f'uJy/duJy'] = self.t['uJy']/self.t[self.dflux_colname]

	def get_median_dflux(self, indices=None):
		if indices is None:
			indices = self.getindices()
		return np.nanmedian(self.t.loc[indices, 'duJy'])

	def get_stdev_flux(self, indices=None):
		self.calcaverage_sigmacutloop('uJy', indices=indices, Nsigma=3.0, median_firstiteration=True)
		return self.statparams['stdev']
	
	def add_noise_to_dflux(self, sigma_extra):
		self.t['duJy_new'] = np.sqrt(self.t['duJy']*self.t['duJy'] + sigma_extra**2)
		self.dflux_colname = 'duJy_new'
		self.calculate_fdf_column()

	def flag_by_control_stats(self, cut:Cut):
		# flag SN measurements according to given bounds
		flag_x2_ix = self.ix_inrange(colnames=['c2_X2norm'], lowlim=cut.params['x2_max'], exclude_lowlim=True)
		flag_stn_ix = self.ix_inrange(colnames=['c2_abs_stn'], lowlim=cut.params['stn_max'], exclude_lowlim=True)
		flag_nclip_ix = self.ix_inrange(colnames=['c2_Nclip'], lowlim=cut.params['Nclip_max'], exclude_lowlim=True)
		flag_ngood_ix = self.ix_inrange(colnames=['c2_Ngood'], uplim=cut.params['Ngood_min'], exclude_uplim=True)
		self.update_mask_column(cut.params['x2_flag'], flag_x2_ix)
		self.update_mask_column(cut.params['stn_flag'], flag_stn_ix)
		self.update_mask_column(cut.params['Nclip_flag'], flag_nclip_ix)
		self.update_mask_column(cut.params['Ngood_flag'], flag_ngood_ix)

		# update mask column with control light curve cut on any measurements flagged according to given bounds
		zero_Nclip_ix = self.ix_equal('c2_Nclip', 0)
		unmasked_ix = self.ix_unmasked('Mask', maskval=cut.params['x2_flag']|cut.params['stn_flag']|cut.params['Nclip_flag']|cut.params['Ngood_flag'])
		self.update_mask_column(cut.params['questionable_flag'], AnotB(unmasked_ix, zero_Nclip_ix))
		self.update_mask_column(cut.flag, AnotB(self.getindices(),unmasked_ix))

	def copy_flags(self, flags_to_copy):
		self.t['Mask'] = self.t['Mask'].astype(np.int32)
		if len(self.t) < 1:
			return
		elif len(self.t) == 1:
			self.t.loc[0,'Mask']= int(self.t.loc[0,'Mask']) | flags_to_copy
		else:
			self.t['Mask'] = np.bitwise_or(self.t['Mask'], flags_to_copy)

	def apply_cut(self, column_name, flag, min_value=None, max_value=None):
		all_ix = self.getindices()
		if min_value:
			kept_ix = self.ix_inrange(colnames=[column_name], lowlim=min_value)
		elif max_value:
			kept_ix = self.ix_inrange(colnames=[column_name], uplim=max_value)
		else:
			raise RuntimeError(f'ERROR: Cannot apply cut without min value ({min_value}) or max value ({max_value}).')
		cut_ix = AnotB(all_ix, kept_ix)
		
		self.update_mask_column(self, flag, cut_ix)

		percent_cut = 100 * len(cut_ix)/len(all_ix)
		return percent_cut

	def update_mask_column(self, flag, indices, remove_old=True):
		if remove_old:
			# remove any old flags of the same value
			self.t['Mask'] = np.bitwise_and(self.t['Mask'].astype(int), ~flag)

		if len(indices) > 1:
			flag_arr = np.full(self.t.loc[indices,'Mask'].shape, flag)
			self.t.loc[indices,'Mask'] = np.bitwise_or(self.t.loc[indices,'Mask'].astype(int), flag_arr)
		elif len(indices) == 1:
			self.t.loc[indices,'Mask'] = int(self.t.loc[indices,'Mask']) | flag

	def drop_extra_columns(self, verbose=False):
		dropcols = []
		for col in ['Noffsetlc', 'uJy/duJy', '__tmp_SN', 'SNR', 'SNRsum', 'SNRsumnorm', 'SNRsim', 'SNRsimsum', 'c2_mean', 'c2_mean_err', 'c2_stdev', 'c2_stdev_err', 'c2_X2norm', 'c2_Ngood', 'c2_Nclip', 'c2_Nmask', 'c2_Nnan', 'c2_abs_stn']:
			if col in self.t.columns:
				dropcols.append(col)
		for col in self.t.columns:
			if re.search('^c\d_',col): 
				dropcols.append(col)

		if len(dropcols) > 0: 
			if verbose:
				print(f'Dropping extra columns ({f"control light curve {str(self.control_index)}" if self.control_index > 0 else "SN light curve"}): ',dropcols)
			self.t.drop(columns=dropcols,inplace=True)

	def check_column_names(self):
		if self.t is None:
			return
		
		for column_name in REQUIRED_COLUMN_NAMES:
			if not column_name in self.t.columns:
				raise RuntimeError(f'ERROR: Missing required column: {column_name}')

	def load(self, input_dir, tnsname):
		filename = get_filename(input_dir, tnsname, self.filt, self.control_index)
		self.load_by_filename(filename)

	def load_by_filename(self, filename):
		self.load_spacesep(filename, delim_whitespace=True, hexcols=['Mask'])
		self.check_column_names()

	def save(self, output_dir, tnsname, indices=None, overwrite=False):
		filename = get_filename(output_dir, tnsname, self.filt, self.control_index)
		self.save_by_filename(filename, indices=indices, overwrite=overwrite)

	def save_by_filename(self, filename, indices=None, overwrite=False):
		self.write(filename=filename, indices=indices, overwrite=overwrite, hexcols=['Mask'])

class AveragedLightCurve(LightCurve):
	def __init__(self, control_index=0, filt='o', mjdbinsize=1.0): 
		LightCurve.__init__(self, control_index, filt) 
		self.mjdbinsize = mjdbinsize

	def load(self, input_dir, tnsname):
		filename = get_filename(input_dir, tnsname, self.filt, self.control_index, self.mjdbinsize)
		self.load_by_filename(filename)
	
	def save(self, output_dir, tnsname, indices=None, overwrite=False):
		filename = get_filename(output_dir, tnsname, self.filt, self.control_index, self.mjdbinsize)
		self.save_by_filename(filename, indices=indices, overwrite=overwrite)