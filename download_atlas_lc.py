#!/usr/bin/env python
'''
Code adapted from Qinan Wang and Armin Rest by Sofia Rest
'''

import configparser, sys, argparse, requests, re, time, json, io, math
import pandas as pd
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord
from collections import OrderedDict
from getpass import getpass
from pdastro import pdastrostatsclass, AorB, AnotB
from atlas_lc import atlas_lc
from plot_atlas_lc import plot_atlas_lc

def RaInDeg(ra):
	s = re.compile('\:')
	if isinstance(ra,str) and s.search(ra):
		A = Angle(ra, u.hour)
	else:
		A = Angle(ra, u.degree)
	return(A.degree)
		   
def DecInDeg(dec):
	A = Angle(dec, u.degree)
	return(A.degree)

class download_atlas_lc:
	def __init__(self):
		# credentials
		self.username = None
		self.password = None
		self.tns_api_key = None

		# input/output
		self.output_dir = None
		self.overwrite = True
		self.flux2mag_sigmalimit = None

		# other
		self.lookbacktime_days = None
		self.controls = False
		self.control_coords = pdastrostatsclass()
		self.radius = None
		self.num_controls = None
		self.closebright_coords = None
		self.closebright_min_dist = None

	# define command line arguments
	def define_args(self, parser=None, usage=None, conflict_handler='resolve'):
		if parser is None:
			parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)
		
		parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')
		
		parser.add_argument('-c','--controls', default=False, action='store_true', help='download control light curves in addition to transient light curve')
		parser.add_argument('-b','--closebright', type=str, default=None, help='comma-separated RA and Dec coordinates of a nearby bright object interfering with the light curve to become center of control light curve circle')
		parser.add_argument('--start_from', type=int, default=1, help='start from a specific control light curve index (useful when program raised error and you don\'t want to have to redownload certain light curves)')')

		"""
		parser.add_argument('-p','--plot', default=False, action='store_true', help='plot light curves and save into PDF file')
		parser.add_argument('--xlim_lower', type=float, default=None, help='if plotting, manually set lower x axis limit to a certain MJD')
		parser.add_argument('--xlim_upper', type=float, default=None, help='if plotting, manually set upper x axis limit to a certain MJD')
		parser.add_argument('--ylim_lower', type=float, default=None, help='if plotting, manually set lower y axis limit to a certain uJy')
		parser.add_argument('--ylim_upper', type=float, default=None, help='if plotting, manually set upper y axis limit to a certain uJy')
		"""

		parser.add_argument('-u','--username', type=str, help='username for ATLAS api')
		parser.add_argument('-a','--tns_api_key', type=str, help='api key to access TNS')
		parser.add_argument('-f','--cfg_filename', default='atlas_lc_settings.ini', type=str, help='file name of ini file with settings for this class')
		parser.add_argument('-l', '--lookbacktime_days', default=None, type=int, help='lookback time in days')
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

		print(f'List of transients to download: {args.tnsnames}')

		self.username = cfg['ATLAS credentials']['username'] if args.username is None else args.username
		print(f'ATLAS username: {self.username}')
		self.password = getpass(prompt='Enter ATLAS password: ')

		self.tns_api_key = cfg['TNS credentials']['api_key'] if args.tns_api_key is None else args.tns_api_key

		self.output_dir = cfg['Input/output settings']['output_dir']
		print(f'Light curve .txt files output directory: {self.output_dir}')
		self.overwrite = not args.dont_overwrite
		print(f'Overwrite existing light curve files: {self.overwrite}')
		self.flux2mag_sigmalimit = int(cfg['Input/output settings']['flux2mag_sigmalimit'])
		print(f'Sigma limit when converting flux to magnitude (magnitudes are limits when dmagnitudes are NaN): {self.flux2mag_sigmalimit}')
		#print(f'Plotting: {args.plot}')
		
		self.lookbacktime_days = args.lookbacktime_days
		if not(self.lookbacktime_days is None): 
			print(f'Lookback time in days: {self.lookbacktime_days}')
		else:
			print('Downloading full light curve(s)')
		
		self.controls = args.controls 
		print(f'Control light curve status: {self.controls}')
		if self.controls:
			self.radius = float(cfg['Control light curve settings']['radius'])
			self.num_controls = int(cfg['Control light curve settings']['num_controls'])
			print(f'# Circle pattern of {self.num_controls} control light curves with radius of {self.radius}\"')
			if not(args.closebright is None):
				self.closebright_coords = [coord.strip() for coord in args.closebright.split(",")]
				if len(self.closebright_coords > 2):
					raise RuntimeError('ERROR: Too many coordinates in --closebright argument!')
				self.closebright_min_dist = float(cfg['Control light curve settings']['closebright_min_dist'])
				print(f'# Close bright object coordinates: RA {self.closebright_coords[0]}, Dec {self.closebright_coords[1]}')
				print(f'# Minimum distance of control light curves from SN location: {self.closebright_min_dist}\"')

	def connect_atlas(self):
		baseurl = 'https://fallingstar-data.com/forcedphot'
		resp = requests.post(url=f"{baseurl}/api-token-auth/",data={'username':self.username,'password':self.password})
		if resp.status_code == 200:
			token = resp.json()['token']
			print(f'Token: {token}')
			headers = {'Authorization':f'Token {token}','Accept':'application/json'}
		else:
			raise RuntimeError(f'ERROR in connect_atlas(): {resp.status_code}')
			print(resp.json())
		return headers

	# API GUIDE: https://fallingstar-data.com/forcedphot/apiguide/
	def get_result(self, ra, dec, headers, lookbacktime_days=None, mjd_max=None):
		if not(lookbacktime_days is None):
			lookbacktime_days = int(Time.now().mjd - lookbacktime_days)
		else:
			lookbacktime_days = int(Time.now().mjd - 1890)
		if not(mjd_max is None):
			mjd_max = int(Time.now().mjd - mjd_max)
		print(f'MJD min: {lookbacktime_days}; MJD max: {mjd_max}')

		baseurl = 'https://fallingstar-data.com/forcedphot'
		task_url = None
		while not task_url:
			with requests.Session() as s:
				resp = s.post(f"{baseurl}/queue/",headers=headers,data={'ra':ra,'dec':dec,'send_email':False,"mjd_min":lookbacktime_days,"mjd_max":mjd_max})
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
					if resp.json()['finishtimestamp']:
						result_url = resp.json()['result_url']
						print(f"Task is complete with results available at {result_url}")
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
			result = s.get(result_url, headers=headers).text
			
		dfresult = pd.read_csv(io.StringIO(result.replace("###", "")), delim_whitespace=True)
		return dfresult

	# get RA and Dec coordinates of control light curves in a circle pattern around SN location and add to control_coords table
	def get_control_coords(self, sn_lc):
		self.control_coords.t = pd.DataFrame(columns=['tnsname','ra','dec','ra_offset','dec_offset','radius','n_detec','n_detec_o','n_detec_c'])

		sn_ra = Angle(RaInDeg(sn_lc.ra), u.degree)
		sn_dec = Angle(DecInDeg(sn_lc.dec), u.degree)

		if self.closebright_coords is None:
			r = Angle(self.radius, u.arcsec)

			# circle pattern center is SN location
			ra_center = Angle(sn_ra.degree, u.degree)
			dec_center = Angle(sn_dec.degree, u.degree)

			# add SN coordinates as first row
			self.control_coords.t = self.control_coords.t.append({'tnsname':sn_lc.tnsname,'ra':sn_ra.degree,'dec':sn_dec.degree,'ra_offset':0,'dec_offset':0,'radius':0,'n_detec':np.nan,'n_detec_o':np.nan,'n_detec_c':np.nan},ignore_index=True)
		else:
			# coordinates of close bright object
			cb_ra = Angle(RaInDeg(closebright_coords[0]), u.degree)
			cb_dec = Angle(DecInDeg(closebright_coords[1]), u.degree)

			# circle pattern radius is distance between SN and bright object
			c1 = SkyCoord(sn_ra, sn_dec, frame='fk5')
			c2 = SkyCoord(cb_ra, cb_dec, frame='fk5')
			sep = c1.separation(c2)
			r = sep.arcsecond

			# circle pattern center is close bright object location
			ra_center = Angle(cb_ra.degree, u.degree)
			dec_center = Angle(cb_dec.degree, u.degree)

			# add SN coordinates as first row; columns like ra_offset, dec_offset, etc. do not apply here
			self.control_coords.t = self.control_coords.t.append({'tnsname':sn_lc.tnsname,'ra':sn_ra.degree,'dec':sn_dec.degree,'ra_offset':np.nan,'dec_offset':np.nan,'radius':np.nan,'n_detec':np.nan,'n_detec_o':np.nan,'n_detec_c':np.nan},ignore_index=True)

		for i in range(self.num_controls):
			angle = Angle(i*360.0/self.num_controls, u.degree)
			
			ra_distance = Angle(r.degree*math.cos(angle.radian), u.degree)
			ra_offset = Angle(ra_distance.degree*(1.0/math.cos(dec_center.radian)), u.degree)
			ra = Angle(ra_center.degree + ra_offset.degree, u.degree)

			dec_offset = Angle(r.degree*math.sin(angle.radian), u.degree)
			dec = Angle(dec_center.degree + dec_offset.degree, u.degree)

			if not(self.closebright_coords is None):
				# check to see if control light curve location is within minimum distance from SN location
				c1 = SkyCoord(sn_ra, sn_dec, frame='fk5')
				c2 = SkyCoord(ra, dec, frame='fk5')
				offset_sep = c1.separation(c2).arcsecond
				if offset_sep < self.closebright_min_dist:
					print(f'Control light curve {i+1:3d} too close to SN location ({offset_sep}\" away) with minimum distance to SN as {self.closebright_min_dist}; skipping control light curve...')
					continue

			# add RA and Dec coordinates to control_coords table
			self.control_coords.t = self.control_coords.t.append({'tnsname':np.nan,'ra':ra.degree,'dec':dec.degree,'ra_offset':ra_offset.degree,'dec_offset':dec_offset.degree,'radius':r,'n_detec':np.nan,'n_detec_o':np.nan,'n_detec_c':np.nan},ignore_index=True)

	# update number of control light curve detections in control_coords table
	def update_control_coords(self, control_lc, control_index):
		o_ix = control_lc.pdastro.ix_equal(colnames=['F'],val='o')
		self.control_coords.t.loc[control_index,'n_detec'] = len(control_lc.pdastro.t)
		self.control_coords.t.loc[control_index,'n_detec_o'] = len(control_lc.pdastro.t.loc[o_ix])
		self.control_coords.t.loc[control_index,'n_detec_c'] = len(control_lc.pdastro.t.loc[AnotB(control_lc.pdastro.getindices(),o_ix)])

	# download a single light curve
	def download_lc(self, args, lc, token):	
		print(f'Downloading forced photometry light curve at {RaInDeg(lc.ra)}, {DecInDeg(lc.dec)} from ATLAS')

		try:
			lc.pdastro.t = self.get_result(RaInDeg(lc.ra), DecInDeg(lc.dec), token, lookbacktime_days=self.lookbacktime_days)
		except Exception as e:
			raise RuntimeError('ERROR in get_result(): '+str(e))

		# sort data by mjd
		lc.pdastro.t = lc.pdastro.t.sort_values(by=['MJD'],ignore_index=True)

		# remove rows with duJy=0 or uJy=Nan
		dflux_zero_ix = lc.pdastro.ix_inrange(colnames='duJy',lowlim=0,uplim=0)
		flux_nan_ix = lc.pdastro.ix_remove_null(colnames='uJy')
		lc.pdastro.t.drop(AorB(dflux_zero_ix,flux_nan_ix))

		lc.pdastro.flux2mag('uJy', 'duJy', 'm', 'dm', zpt=23.9, upperlim_Nsigma=self.flux2mag_sigmalimit)

		return lc

	# download SN light curve and, if necessary, control light curves, then save
	def download_lcs(self, args, tnsname, token):
		lc = atlas_lc(tnsname=tnsname)
		print(f'\nCOMMENCING LOOP FOR SN {lc.tnsname}\n')
		lc.get_tns_data(self.tns_api_key)
		lc = self.download_lc(args, lc, token)
		lc.save_lc(self.output_dir, overwrite=self.overwrite)

		"""
		if args.plot:
			plot = plot_atlas_lc(tnsname=lc.tnsname, output_dir=self.output_dir, args=args)
			plot.set(lc=lc, filt=filt)
			plot.plot_og_lc(separate_baseline=False)
		"""

		if args.controls:
			print('Control light curve downloading set to True')
			
			self.get_control_coords(lc)
			print('Control light curve coordinates calculated: \n',self.control_coords.t[['tnsname','ra','dec','ra_offset','dec_offset','radius']])

			# download control light curves
			for control_index in range(args.start_from,len(self.control_coords.t)):
				print(f'\nControl light curve {control_index:03d}:\n',self.control_coords.t.loc[control_index,['ra','dec','ra_offset','dec_offset','radius']].T)
				control_lc = atlas_lc(ra=self.control_coords.t.loc[control_index,'ra'], dec=self.control_coords.t.loc[control_index,'dec'])
				control_lc = self.download_lc(args, control_lc, token)
				self.update_control_coords(control_lc, control_index)
				lc.add_control_lc(control_lc.pdastro)
				lc.save_control_lc(self.output_dir, control_index, overwrite=self.overwrite)

			# save control_coords table
			self.control_coords.write(filename=f'{self.output_dir}/{lc.tnsname}_control_coords.txt', overwrite=self.overwrite)

			"""
			if args.plot:
				plot.plot_og_control_lcs()
			"""
		
		#lc.save(self.output_dir, overwrite=self.overwrite)

	# loop through each SN given and download light curves
	def download_loop(self):
		args = self.define_args().parse_args()
		self.load_settings(args)

		print('\nConnecting to ATLAS API...')
		token = self.connect_atlas()
		if token is None: 
			raise RuntimeError('ERROR in connect_atlas(): No token header!')

		for obj_index in range(0,len(args.tnsnames)):
			self.download_lcs(args, args.tnsnames[obj_index], token)

if __name__ == "__main__":
	download_atlas_lc = download_atlas_lc()
	download_atlas_lc.download_loop()
