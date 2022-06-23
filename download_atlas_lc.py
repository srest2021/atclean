import configparser, sys, argparse, requests, re, time
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord
from pdastro import pdastrostatsclass, AorB

class lc:
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
			raise RuntimeError('ERROR: \n'+str(e))

		self.ra = json_data['data']['reply']['ra']
		self.dec = json_data['data']['reply']['dec']

		discoverydate = json_data['data']['reply']['discoverydate']
		date = list(discoverydate.partition(' '))[0]
		time = list(discoverydate.partition(' '))[2]
		dateobjects = Time(date+"T"+time, format='isot', scale='utc')
		self.discdate = dateobjects.mjd

	# get a light curve filename for saving
	def get_filename(self, filt, control_index):
		if not self.is_averaged:
			filename = f'{self.output_dir}/{self.tnsname}_i{control_index:3d}.{filt}.lc.txt'
		else:
			filename = f'{self.output_dir}/{self.tnsname}_i{control_index:3d}.{filt}.{self.mjdbinsize:0.2f}days.lc.txt'
		print(f'# Filename: {filename}')
		return filename

	# save SN light curve and, if necessary, control light curves
	def save(self, overwrite=True):
		print('Saving SN light curve')
		o_ix = self.pdastro.ix_equal(colnames=['F'],val='o')
		self.pdastro.write(filename=get_filename('o',i), indices=o_ix, overwrite=overwrite)
		self.pdastro.write(filename=get_filename('c',i), indices=AnotB(self.pdastro.getindices(),o_ix), overwrite=overwrite)

		if len(self.lcs) > 0:
			print('Saving control light curves')
			for i in range(1,len(self.lcs)):
				for filt in ['c','o']:
					filt_ix = self.lcs[i].pdastro.ix_equal(colnames=['F'],val=filt)
					self.lcs[i].pdastro.write(filename=get_filename(filt,i), indices=filt_ix, overwrite=overwrite)

	# add already downloaded control light curve to control light curve dictionary
	def add_control_lc(self, control_lc):
		self.lcs[len(self.lcs)+1] = control_lc

class download_atlas_lc:
	def __init__(self):
		# credentials
		self.username = None
		self.password = None
		self.tns_api_key = None

		# output
		self.output_dir = None
		self.overwrite = True

		# other
		self.lookbacktime_days = None
		self.controls = False
		self.control_coords = pdastrostatsclass()
		self.radius = None
		self.num_controls = None

	# define command line arguments
	def define_args(self, parser=None, usage=None, conflict_handler='resolve'):
		if parser is None:
			parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)
		
		parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')
		parser.add_argument('-c','--controls', default=False, action='store_true', help='download control light curves in addition to transient light curve')
		parser.add_argument('-u','--username', type=str, help='username for ATLAS api')
		parser.add_argument('-p','--password', type=str, default=None, help='password for ATLAS api')
		parser.add_argument('-a','--tns_api_key', type=str, help='api key to access TNS')
		parser.add_argument('-f','--cfg_filename', default='atlaslc.ini', type=str, help='file name of ini file with settings for this class')
		parser.add_argument('-l', '--lookbacktime_days', default=None, type=int, help='lookback time in days')
		parser.add_argument('-o', '--dont_overwrite', default=False, action='store_true', help='overwrite existing file with same file name')
		
		return parser

	# load config settings from file and reconcile with command arguments
	def load_settings(self, args):
		cfg = configparser.ConfigParser()
		cfg.read(args.cfg_filename)

		self.username = cfg['ATLAS credentials']['username'] if args.username is None else args.username
		if args.password is None:
			raise RuntimeError('ERROR: please provide ATLAS password using --password argument!')
		self.password = args.password
		self.tns_api_key = cfg['TNS credentials']['api_key'] if args.tns_api_key is None else args.tns_api_key

		self.output_dir = cfg['Input/output settings']['output_dir']

		self.overwrite = not args.dont_overwrite
		self.lookbacktime_days = args.lookbacktime_days
		self.controls = args.controls 
		if self.controls:
			self.radius = cfg['Control light curve settings']['radius']
			self.num_controls = cfg['Control light curve settings']['num_controls']

	def connect_atlas(self):
		resp = requests.post(url=f"{self.baseurl}/api-token-auth/",data={'username':self.username,'password':self.password})
		if resp.status_code == 200:
			token = resp.json()['token']
			print(f'Token: {token}')
			headers = {'Authorization':f'Token {token}','Accept':'application/json'}
		else:
			print(f'ERROR: {resp.status_code}')
			print(resp.json())
		return headers

	# API GUIDE: https://fallingstar-data.com/forcedphot/apiguide/
	def get_result(self, ra, dec, headers, lookbacktime_days=None, mjd_max=None):
		if not(lookbacktime_days is None): lookbacktime_days = int(Time.now().mjd - lookbacktime_days)
		else: lookbacktime_days = int(Time.now().mjd - 1890)
		if not(mjd_max is None): mjd_max = int(Time.now().mjd - mjd_max)
		print('MJD min: ',lookbacktime_days,'. MJD max: ',mjd_max)

		task_url = None
		while not task_url:
			with requests.Session() as s:
				resp = s.post(f"{self.baseurl}/queue/",headers=headers,data={'ra':ra,'dec':dec,'send_email':False,"mjd_min":lookbacktime_days,"mjd_max":mjd_max})
				if resp.status_code == 201:  # successfully queued
					task_url = resp.json()['url']
					print(f'The task URL is {task_url}')
					print(f"Waiting for job to start (queued at {resp.json()['timestamp']})")
				elif resp.status_code == 429:  # throttled
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
		while not result_url:
			with requests.Session() as s:
				resp = s.get(task_url, headers=headers)
				if resp.status_code == 200:  # HTTP OK
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

		# add SN coordinates as first row
		self.control_coords.t = self.control_coords.t.append({'tnsname':sn_lc.tnsname,'ra':sn_lc.ra,'dec':sn_lc.dec,'ra_offset':0,'dec_offset':0,'radius':0,'n_detec':np.nan,'n_detec_o':np.nan,'n_detec_c':np.nan},ignore_index=True)

		r = Angle(self.radius, u.arcsec)
		ra_center = Angle(sn_lc.ra.degree, u.degree)
		dec_center = Angle(sn_lc.dec.degree, u.degree)

		for i in range(self.num_controls):
			angle = Angle(i*360.0/self.num_controls, u.degree)
			
			ra_distance = Angle(r.degree*math.cos(angle.radian), u.degree)
			ra_offset = Angle(ra_distance.degree*(1.0/math.cos(dec_center.radian)), u.degree)
			ra = Angle(ra_cecter.degree + ra_offset.degree, u.degree)

			dec_offset = Angle(r.degree*math.sin(angle.radian), u.degree)
			dec = Angle(dec_center.degree + dec_offset.degree, u.degree)

			# add RA and Dec coordinates to control_coords table
			self.control_coords.t = self.control_coords.t.append({'tnsname':np.nan,'ra':ra,'dec':dec,'ra_offset':ra_offset,'dec_offset':dec_offset,'radius':r,'n_detec':np.nan,'n_detec_o':np.nan,'n_detec_c':np.nan},ignore_index=True)

	# update number of control light curve detections in control_coords table
	def update_control_coords(self, control_lc, control_index):
		o_ix = control_lc.pdastro.ix_equal(colnames=['F'],val='o')
		self.control_coords.t.loc[control_index,'n_detec'] = len(control_lc.pdastro.t)
		self.control_coords.t.loc[control_index,'n_detec_o'] = len(control_lc.pdastro.t.loc[o_ix])
		self.control_coords.t.loc[control_index,'n_detec_c'] = len(control_lc.pdastro.t.loc[AnotB(control_lc.pdastro.getindices(),o_ix)])

	# download a single light curve
	def download_lc(self, lc):	
		print(f'Downloading forced photometry light curve at {lc.ra}, {lc.dec} from ATLAS')

		lc.pdastro.t = get_result(lc.ra, lc.dec, token, lookbacktime_days=args.lookbacktime_days)

		# sort data by mjd
		lc.pdastro.t = lc.pdastro.t.sort_values(by=['MJD'],ignore_index=True)

		# remove rows with duJy=0 or uJy=Nan
		dflux_zero_ix = lc.pdastro.ix_inrange(colnames='duJy',lowlim=0,uplim=0)
		flux_nan_ix = lc.pdastro.ix_remove_null(colnames='uJy')
		lc.pdastro.t.drop(AorB(dflux_zero_ix,flux_nan_ix))

		return lc

	# download SN light curve and, if necessary, control light curves, then save
	def download_lcs(self, args, obj_index, token):
		lc = lc(tnsname=args.tnsnames[obj_index])
		print(f'Commencing download loop for SN {lc.tnsname}')
		lc.get_tns_data(self.tns_api_key)
		lc = download_lc(lc)

		if args.controls:
			print('Control light curve downloading set to True')
			
			get_control_coords(lc)
			print('Control light curve coordinates calculated: \n',self.control_coords.t[['tnsname','ra','dec','ra_offset','dec_offset','radius']])

			# download control light curves
			for i in range(1,len(self.control_coords.t)):
				control_lc = lc(ra=self.control_coords.t.loc[i,'ra'], dec=self.control_coords.t.loc[i,'dec'])
				control_lc = download_lc(control_lc)
				update_control_coords(control_lc, i)
				lc.add_control_lc(control_lc)

			# save control_coords table
			self.control_coords.write(filename=f'{self.output_dir}/{lc.tnsname}_control_coords.txt', overwrite=self.overwrite)
		
		lc.save()

	# loop through each SN given and download light curves
	def download_loop(self):
		args = self.define_args().parse_args()

		load_settings(args)

		print('Connecting to ATLAS API...')
		token = connect_atlas()
		if token is None: 
			raise RuntimeError('ERROR: No token header!')

		for obj_index in range(0,len(args.tnsnames)):
			download_lcs(args, obj_index, token)

