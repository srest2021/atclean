#!/usr/bin/env python
'''
Code adapted from Qinan Wang and Armin Rest by Sofia Rest
'''

import configparser, sys, argparse, requests, re, time, io, math, os
import pandas as pd
import numpy as np
from getpass import getpass
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord

from pdastro import pdastrostatsclass, AorB, AnotB
from atlas_lc import atlas_lc

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
		self.tns_id = None
		self.bot_name = None

		# input/output
		self.output_dir = None
		self.snlist_filename = None
		self.snlist = None
		self.overwrite = True
		self.flux2mag_sigmalimit = None

		self.mjd_min = None
		self.mjd_max = None

		# other
		self.controls = False
		self.control_coords = pdastrostatsclass()
		self.radius = None
		self.num_controls = None
		self.closebright = False
		self.closebright_coords = None
		self.closebright_min_dist = None
		self.control_coords_filename = None

	# define command line arguments
	def define_args(self, parser=None, usage=None, conflict_handler='resolve'):
		if parser is None:
			parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)
		
		parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')
		parser.add_argument('--coords', type=str, default=None, help='comma-separated RA and Dec of SN light curve to download')
		parser.add_argument('--discdate', type=str, default=None, help='SN discovery date in MJD')

		parser.add_argument('-c','--controls', default=False, action='store_true', help='download control light curves in addition to transient light curve')
		parser.add_argument('--closebright', type=str, default=None, help='comma-separated RA and Dec coordinates of a nearby bright object interfering with the light curve to become center of control light curve circle')
		parser.add_argument('--ctrl_coords', type=str, default=None, help='file name of text file in output_dir containing table of control light curve coordinates')

		parser.add_argument('-u','--username', type=str, help='username for ATLAS api')
		parser.add_argument('-a','--tns_api_key', type=str, help='api key to access TNS')
		parser.add_argument('-f','--cfg_filename', default='params.ini', type=str, help='file name of ini file with settings for this class')
		parser.add_argument('-l', '--lookbacktime_days', default=None, type=int, help='lookback time (days)')
		parser.add_argument('--mjd_min', default=None, type=float, help='minimum MJD to download')
		parser.add_argument('--mjd_max', default=None, type=float, help='maximum MJD to download')
		parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite existing file with same file name')

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

		# ATLAS credentials
		self.username = cfg['ATLAS credentials']['username'] if args.username is None else args.username
		print(f'ATLAS username: {self.username}')
		self.password = getpass(prompt='Enter ATLAS password: ')

		# TNS credentials
		self.tns_api_key = cfg['TNS credentials']['api_key'] if args.tns_api_key is None else args.tns_api_key
		self.tns_id = cfg['TNS credentials']['tns_id']
		self.bot_name = cfg['TNS credentials']['bot_name']

		# data output directory
		self.output_dir = cfg['Input/output settings']['output_dir']
		print(f'Light curve .txt files output directory: {self.output_dir}')

		# attempt to load snlist.txt; if does not exist, create new snlist table
		self.snlist_filename = f'{self.output_dir}/{cfg["Input/output settings"]["snlist_filename"]}'
		if os.path.exists(self.snlist_filename):
			self.snlist = pdastrostatsclass()
			print(f'Loading SN list at {self.snlist_filename}')
			self.snlist.load_spacesep(self.snlist_filename, delim_whitespace=True)
		else:
			self.snlist = pdastrostatsclass(columns=['tnsname', 'ra', 'dec', 'discovery_date', 'closebright_ra', 'closebright_dec'])

		# overwrite existing files?
		self.overwrite = bool(args.overwrite)
		print(f'Overwrite existing light curve files: {self.overwrite}')

		# flux2mag sigma limit
		self.flux2mag_sigmalimit = int(cfg['Input/output settings']['flux2mag_sigmalimit'])
		print(f'Sigma limit when converting flux to magnitude (magnitudes are limits when dmagnitudes are NaN): {self.flux2mag_sigmalimit}')
		
		# mjd min
		if not(args.lookbacktime_days is None) and not(args.mjd_min is None):
			raise RuntimeError('ERROR: Cannot specify minimum MJD and lookback time at the same time! Choose one argument only')
		elif not(args.lookbacktime_days is None):
			print(f'Lookback time (days): {args.lookbacktime_days}')
			self.mjd_min = float(Time.now().mjd - args.lookbacktime_days)
		elif not(args.mjd_min is None):
			self.mjd_min = float(args.mjd_min)
		else:
			self.mjd_min = 50000.0
		print(f'Min MJD: {self.mjd_min:0.2f}')
		# mjd max
		if not(args.mjd_max is None): 
			self.mjd_max = float(args.mjd_max)
		else:
			self.mjd_max = float(Time.now().mjd)
		print(f'Max MJD: {self.mjd_max:0.2f}')
		
		# control lcs
		self.controls = bool(args.controls)
		print(f'Control light curve status: {self.controls}')
		if self.controls:
			if args.ctrl_coords is None:
				self.radius = float(cfg['Control light curve settings']['radius'])
				self.num_controls = int(cfg['Control light curve settings']['num_controls'])
				print(f'# Circle pattern of {self.num_controls} control light curves with radius of {self.radius}\"')
				
				if args.closebright:
					self.closebright = True
					if len(args.tnsnames) > 1:
						raise RuntimeError('ERROR: Only one SN allowed when specifying close bright object coordinates in command line!')
					
					self.closebright_coords = [coord.strip() for coord in args.closebright.split(",")]
					if len(self.closebright_coords) > 2:
						raise RuntimeError('ERROR: Too many coordinates in --closebright argument!')
					print(f'# Close bright object coordinates: RA {self.closebright_coords[0]}, Dec {self.closebright_coords[1]}')
					
					self.closebright_min_dist = float(cfg['Control light curve settings']['closebright_min_dist'])
					print(f'# Minimum distance of control light curves from SN location: {self.closebright_min_dist}\"')
			else:
				self.control_coords_filename = f'{self.output_dir}/{args.ctrl_coords}'
				print(f'# Path name of txt file with control light curve coordinates: {self.control_coords_filename}')
			
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
	def get_result(self, ra, dec, headers, mjd_min, mjd_max):
		if mjd_min > mjd_max:
			raise RuntimeError(f'ERROR: max MJD {mjd_max} less than min MJD {mjd_min}')
			sys.exit()
		else:
			print(f'Min MJD: {mjd_min}; max MJD: {mjd_max}')
		
		baseurl = 'https://fallingstar-data.com/forcedphot'
		task_url = None
		while not task_url:
			with requests.Session() as s:
				resp = s.post(f"{baseurl}/queue/",headers=headers,data={'ra':ra,'dec':dec,'send_email':False,"mjd_min":mjd_min,"mjd_max":mjd_max})
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
				print('WARNING: Empty light curve--no data within this MJD range...')
				dfresult = pd.DataFrame(columns=['MJD','m','dm','uJy','duJy','F','err','chi/N','RA','Dec','x','y','maj','min','phi','apfit','Sky','ZP','Obs','Mask'])
			else:
				result = s.get(result_url, headers=headers).text
				dfresult = pd.read_csv(io.StringIO(result.replace("###", "")), delim_whitespace=True)
		
		return dfresult

	def get_filt_lens(self, sn_lc):
		if len(sn_lc.lcs) > 0: # SN lc has been downloaded
			return sn_lc.get_filt_lens()
		else: # SN lc has not been downloaded, so temporarily load both c and o lcs to get indices
			temp_o = pdastrostatsclass()
			temp_o.load_spacesep(sn_lc.get_filename('o', 0, self.output_dir), delim_whitespace=True)
			o_len = len(temp_o.t)

			temp_c = pdastrostatsclass()
			temp_c.load_spacesep(sn_lc.get_filename('c', 0, self.output_dir), delim_whitespace=True)
			c_len = len(temp_c.t)
			return o_len, c_len

	# convert ra string to angle
	def ra_str2angle(self, ra):
		return Angle(RaInDeg(ra), u.degree)

	# convert dec string to angle
	def dec_str2angle(self, dec):
		return Angle(DecInDeg(dec), u.degree)

	# get distance between 2 locations specified by ra and dec angles
	def get_distance(self, ra1, dec1, ra2, dec2):
		c1 = SkyCoord(ra1, dec1, frame='fk5')
		c2 = SkyCoord(ra2, dec2, frame='fk5')
		return c1.separation(c2)

	def read_coords_from_file(self, sn_ra, sn_dec):
		print(f'Loading control coordinates file at {self.control_coords_filename}...')
		temp = pdastrostatsclass()
		temp.load_spacesep(self.control_coords_filename, delim_whitespace=True)

		if len(temp.t) < 1:
			raise RuntimeError(f'ERROR: no rows in {self.control_coords_filename}--must have at least one control lc location!')
		self.num_controls = len(temp.t)

		for i in range(len(temp.t)):
			ra = self.ra_str2angle(temp.t.loc[i,'ra'])
			dec = self.dec_str2angle(temp.t.loc[i,'dec'])
			r = self.get_distance(sn_ra, sn_dec, ra, dec)

			self.control_coords.newrow({'tnsname':np.nan,
										'control_id':i,
										'ra': f'{ra.degree:0.8f}',
										'dec': f'{dec.degree:0.8f}',
										'ra_offset': np.nan,
										'dec_offset': np.nan,
										'radius':r.arcsecond,
										'n_detec':np.nan,
										'n_detec_o':np.nan,
										'n_detec_c':np.nan})

		with pd.option_context('display.float_format', '{:,.8f}'.format):
			print('Control light curve coordinates read from txt file: \n',self.control_coords.t[['tnsname','control_id','ra','dec','ra_offset','dec_offset','radius']])

	# get RA and Dec coordinates of control light curves in a circle pattern around SN location and add to control_coords table
	def get_control_coords(self, sn_lc):
		self.control_coords.t = pd.DataFrame(columns=['tnsname','control_id','ra','dec','ra_offset','dec_offset','radius','n_detec','n_detec_o','n_detec_c'])

		sn_ra = self.ra_str2angle(sn_lc.ra)
		sn_dec = self.dec_str2angle(sn_lc.dec)
		o_len, c_len = self.get_filt_lens(sn_lc)

		# set first row of control_coords table according to closebright status
		if not self.closebright: # pattern around SN location
			if self.control_coords_filename is None:
				r = Angle(self.radius, u.arcsec)

			# circle pattern center is SN location
			ra_center = sn_ra
			dec_center = sn_dec

			# add SN coordinates as first row
			self.control_coords.newrow({'tnsname':sn_lc.tnsname,
										'control_id':0,
										'ra': f'{sn_ra.degree:0.8f}',
										'dec': f'{sn_dec.degree:0.8f}',
										'ra_offset':0,
										'dec_offset':0,
										'radius':0,
										'n_detec':o_len+c_len,
										'n_detec_o':o_len,
										'n_detec_c':c_len})
		
		else: # pattern around close bright object
			# coordinates of close bright object
			cb_ra = self.ra_str2angle(self.closebright_coords[0])
			cb_dec = self.dec_str2angle(self.closebright_coords[1])

			# circle pattern radius is distance between SN and bright object
			r = self.get_distance(sn_ra, sn_dec, cb_ra, cb_dec).arcsecond

			# circle pattern center is close bright object location
			ra_center = cb_ra
			dec_center = cb_dec

			# add SN coordinates as first row; columns like ra_offset, dec_offset, etc. do not apply here
			self.control_coords.newrow({'tnsname':sn_lc.tnsname,
										'control_id':0,
										'ra': f'{sn_ra.degree:0.8f}',
										'dec': f'{sn_dec.degree:0.8f}',
										'ra_offset':np.nan,
										'dec_offset':np.nan,
										'radius':np.nan,
										'n_detec':o_len+c_len,
										'n_detec_o':o_len,
										'n_detec_c':c_len},ignore_index=True)

		if not(self.control_coords_filename is None): # read in control lc coordinates from txt file
			self.read_coords_from_file(sn_ra, sn_dec)
		
		else: # calculate control lc coordinates
			for i in range(self.num_controls):
				angle = Angle(i*360.0 / self.num_controls, u.degree)
				
				ra_distance = Angle(r.degree * math.cos(angle.radian), u.degree)
				ra_offset = Angle(ra_distance.degree * (1.0/math.cos(dec_center.radian)), u.degree)
				ra = Angle(ra_center.degree + ra_offset.degree, u.degree)

				dec_offset = Angle(r.degree * math.sin(angle.radian), u.degree)
				dec = Angle(dec_center.degree + dec_offset.degree, u.degree)

				if self.closebright: # check to see if control light curve location is within minimum distance from SN location
					offset_sep = self.get_distance(sn_ra, sn_dec, ra, dec).arcsecond
					if offset_sep < self.closebright_min_dist:
						print(f'Control light curve {i+1:3d} too close to SN location ({offset_sep}\" away) with minimum distance to SN as {self.closebright_min_dist}; skipping control light curve...')
						continue

				# add RA and Dec coordinates to control_coords table
				self.control_coords.newrow({'tnsname':np.nan,
											'control_id':i,
											'ra': f'{ra.degree:0.8f}',
											'dec': f'{dec.degree:0.8f}',
											'ra_offset': f'{ra_offset.degree:0.8f}',
											'dec_offset': f'{dec_offset.degree:0.8f}',
											'radius':r,
											'n_detec':np.nan,
											'n_detec_o':np.nan,
											'n_detec_c':np.nan})

			with pd.option_context('display.float_format', '{:,.8f}'.format):
				print('Control light curve coordinates calculated: \n',self.control_coords.t[['tnsname','control_id','ra','dec','ra_offset','dec_offset','radius']])

	# update number of control light curve detections in control_coords table
	def update_control_coords(self, lc, control_index):
		o_ix = lc.lcs[control_index].ix_equal(colnames=['F'],val='o')
		self.control_coords.t.loc[control_index,'n_detec'] = len(lc.lcs[control_index].t)
		self.control_coords.t.loc[control_index,'n_detec_o'] = len(o_ix)
		self.control_coords.t.loc[control_index,'n_detec_c'] = len(AnotB(lc.lcs[control_index].getindices(),o_ix))

	def get_lc_data(self, args, lc, snlist_index):
		# get RA, Dec, and discovery date from command line if available
		if not(args.coords is None):
			print('Getting RA and Dec coordinates from command line argument --coords...')
			if len(args.tnsnames) > 1:
				raise RuntimeError('ERROR: Only one SN allowed when specifying SN coordinates in command line!')
			coords = [coord.strip() for coord in args.coords.split(",")]
			if len(coords) > 2:
				raise RuntimeError('ERROR: Too many coordinates in --coords argument! Please provide comma-separated RA and Dec.')
			if len(coords) < 2:
				raise RuntimeError('ERROR: Too few coordinates in --coords argument! Please provide comma-separated RA and Dec.')
			lc.ra = coords[0]
			lc.dec = coords[1]
		if not(args.discdate is None):
			print('Getting discovery date from command line argument --discdate...')
			if len(args.tnsnames) > 1:
				raise RuntimeError('ERROR: Only one SN allowed when specifying SN discovery date in command line!')
			lc.discdate = args.discdate

		# if RA, Dec, or discovery date not provided in command line, check in snlist.txt or try using TNS API
		if args.coords is None or args.discdate is None:
			# check if no existing row in snlist.txt; else no need to add new row
			if snlist_index == -1:
				# check if we can query TNS
				if self.tns_api_key == 'None':
					# all options exhausted
					raise RuntimeError('ERROR: No TNS API key provided, no corresponding SN entry in snlist.txt, and not enough information provided in arguments! Please provide RA, Dec, and discovery date in --coords and --discdate arguments.')
				else:
					# get RA, Dec, and/or discovery date from TN 
					lc.get_tns_data(self.tns_api_key, self.tns_id, self.bot_name)
			else:
				print(f'Getting RA, Dec, and/or discovery date from existing row in SN list...')

		# add RA, Dec, and discovery date to new row in snlist.txt if no existing row for this SN
		if snlist_index == -1:
			# add row to snlist.txt
			self.snlist.newrow({'tnsname':lc.tnsname, 'ra':lc.ra, 'dec':lc.dec, 'discovery_date':lc.discdate, 'closebright_ra':np.nan, 'closebright_dec':np.nan})
			snlist_index = len(self.snlist.t)-1

		# fill in missing information from snlist.txt
		if lc.ra is None: lc.ra = self.snlist.t.loc[snlist_index,'ra']
		if lc.dec is None: lc.dec = self.snlist.t.loc[snlist_index,'dec']
		if lc.discdate is None: lc.discdate = self.snlist.t.loc[snlist_index,'discovery_date']

		if self.closebright and self.closebright_coords is None:
			print(f'Getting close bright object coordinates from SN list at {self.snlist_filename}...')
			if not(np.isnan(self.snlist.t.loc[snlist_index,'closebright_ra'])) and not(np.isnan(self.snlist.t.loc[snlist_index,'closebright_dec'])):
				self.closebright_coords[0] = self.snlist.t.loc[snlist_index,'closebright_ra']
				self.closebright_coords[1] = self.snlist.t.loc[snlist_index,'closebright_dec']
			else:
				raise RuntimeError(f'ERROR: Close bright object coordinates given in SN list file at {self.snlist_filename} are not valid!')
		
		output = f'RA: {lc.ra}, Dec: {lc.dec}, discovery date: {lc.discdate}'
		if self.closebright: 
			output += f', closebright RA: {self.closebright_coords[0]}, closebright Dec: {self.closebright_coords[1]}'
		print(output)

		return lc

	# download a single light curve
	def download_lc(self, args, token, lc, ra, dec, control_index=0):	
		print(f'Downloading forced photometry light curve at {RaInDeg(ra):0.8f}, {DecInDeg(dec):0.8f} from ATLAS')
		lc.lcs[control_index] = pdastrostatsclass()

		while(True):
			try:
				lc.lcs[control_index].t = self.get_result(RaInDeg(ra), DecInDeg(dec), token, self.mjd_min, self.mjd_max)
				break
			except Exception as e:
				print('Exception caught: '+str(e))
				print('Trying again in 20 seconds! Waiting...')
				time.sleep(20)
				continue

		# sort data by mjd
		lc.lcs[control_index].t = lc.lcs[control_index].t.sort_values(by=['MJD'],ignore_index=True)

		# remove rows with duJy=0 or uJy=Nan
		dflux_zero_ix = lc.lcs[control_index].ix_inrange(colnames='duJy',lowlim=0,uplim=0)
		flux_nan_ix =lc.lcs[control_index].ix_is_null(colnames='uJy')
		print('\nDeleting %d rows with "duJy"==0 or "uJy"==NaN...' % (len(dflux_zero_ix) + len(flux_nan_ix)))
		if len(AorB(dflux_zero_ix,flux_nan_ix)) > 0:
			lc.lcs[control_index].t = lc.lcs[control_index].t.drop(AorB(dflux_zero_ix,flux_nan_ix))
			 
		lc.lcs[control_index].flux2mag('uJy', 'duJy', 'm', 'dm', zpt=23.9, upperlim_Nsigma=self.flux2mag_sigmalimit)

		return lc

	# download SN light curve and, if necessary, control light curves, then save
	def download_lcs(self, args, tnsname, token, snlist_index):
		lc = atlas_lc(tnsname=tnsname)
		print(f'\nCOMMENCING LOOP FOR SN {lc.tnsname}\n')
		
		lc = self.get_lc_data(args, lc, snlist_index)

		# only download light curve if overwriting existing files
		if not(self.overwrite) and lc.exists(self.output_dir, 'o') and lc.exists(self.output_dir, 'c'):
			print(f'SN light curve files already exist and overwrite is set to {self.overwrite}! Skipping download...')
		else:
			lc = self.download_lc(args, token, lc, lc.ra, lc.dec)
			lc.save_lc(self.output_dir, overwrite=self.overwrite)

		if args.controls:
			print('Control light curve downloading set to True')
			
			self.get_control_coords(lc)

			# download control light curves
			for control_index in range(1,len(self.control_coords.t)):
				# only download control light curve if overwriting existing files
				if not(self.overwrite) and lc.exists(self.output_dir, 'o', control_index=control_index) and lc.exists(self.output_dir, 'c', control_index=control_index):
					print(f'Control light curve {control_index:03d} files already exist and overwrite is set to {self.overwrite}! Skipping download...')
				else:
					print(f'\nDownloading control light curve {control_index:03d}...')
					lc = self.download_lc(args, token, lc, 
										  ra=self.control_coords.t.loc[control_index,'ra'], 
										  dec=self.control_coords.t.loc[control_index,'dec'], 
										  control_index=control_index)
					self.update_control_coords(lc, control_index)
					lc.save_lc(self.output_dir, control_index=control_index, overwrite=self.overwrite)

			# save control_coords table
			#self.control_coords.write()
			#self.control_coords.default_formatters = {'ra':'{:.8f}'.format,'dec':'{:.8d}'.format,'ra_offset':'{:.8d}'.format,'dec_offset':'{:.8d}'.format}
			self.control_coords.write(filename=f'{self.output_dir}/{lc.tnsname}/controls/{lc.tnsname}_control_coords.txt', overwrite=self.overwrite)

	# loop through each SN given and download light curves
	def download_loop(self):
		args = self.define_args().parse_args()
		self.load_settings(args)

		print('\nConnecting to ATLAS API...')
		token = self.connect_atlas()
		if token is None: 
			raise RuntimeError('ERROR in connect_atlas(): No token header!')

		for obj_index in range(len(args.tnsnames)):
			snlist_index = -1
			snlist_ix = self.snlist.ix_equal(colnames=['tnsname'],val=args.tnsnames[obj_index])
			# check if SN information exists in snlist.txt
			if len(snlist_ix) > 0:
				if len(snlist_ix > 1):
					# drop duplicate rows
					self.snlist.t.drop(snlist_ix[1:])
				snlist_index = snlist_ix[0]

			self.download_lcs(args, args.tnsnames[obj_index], token, snlist_index)

		# save snlist.txt with any new rows
		print(f'\nSaving SN list at {self.snlist_filename}')
		self.snlist.write(self.snlist_filename)

if __name__ == "__main__":
	download_atlas_lc = download_atlas_lc()
	download_atlas_lc.download_loop()
