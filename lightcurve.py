#!/usr/bin/env python

from typing import Dict, Type
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
  def __init__(self, string:str|None=None):
    self.string:str|None = string
    self.degrees:u.degree|None = self.to_degrees()

  def to_degrees(self):
    if self.string is None:
      return None
    
    s = re.compile('\:')
    if isinstance(self.string,str) and s.search(self.string):
      A = Angle(self.string, u.hour)
    else:
      A = Angle(self.string, u.degree)
    return A.degree
  
  # def __str__(self):
  #   return f'{self.degrees}'
  
class Dec:
  def __init__(self, string:str|None=None):
    self.string:str|None = string
    self.degrees:u.degree|None = self.to_degrees()

  def to_degrees(self):
    if self.string is None:
      return None
    
    A = Angle(self.string, u.degree)
    return A.degree
  
  # def __str__(self):
  #   return f'{self.degrees}'
  
class Coordinates:
  def __init__(self, ra:str|None=None, dec:str|None=None):
    self.ra:RA = RA(ra)
    self.dec:Dec = Dec(dec)

  def set_RA(self, ra):
    self.ra = RA(ra)
  
  def set_Dec(self, dec):
    self.dec = Dec(dec)

  def is_empty(self) -> bool:
    return self.ra.string is None or self.dec.string is None

  def __str__(self):
    if self.is_empty():
      raise RuntimeError(f'ERROR: Coordinates are empty and cannot be printed.')
    return f'RA {self.ra.degrees:0.8f}, Dec {self.dec.degrees:0.8f}'
  
def get_filename(output_dir, tnsname, filt='o', control_index=0, mjdbinsize=None):
  filename = f'{output_dir}/{tnsname}'
  
  if control_index != 0:
    filename += '/controls'
  
  filename += f'/{tnsname}'
  
  if control_index != 0:
    filename += f'_i{control_index:03d}'
  
  filename += f'.{filt}'
  
  if mjdbinsize:
    filename += f'.{mjdbinsize:0.2f}days'
    
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
class SnInfo():
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
      self.t = pd.DataFrame(columns=['tnsname', 'ra', 'dec', 'mjd0', 'closebright_ra', 'closebright_dec'])

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
    self.t.to_string(self.filename)

  def __str__(self):
    print(self.t.to_string())
  
"""
LIGHT CURVES
"""

class Supernova:
  def __init__(self, tnsname:str=None, ra:str=None, dec:str=None, mjd0:float|None=None):
    self.tnsname = tnsname
    self.num_controls = 0
    self.coords:Coordinates = Coordinates(ra,dec)
    self.mjd0 = mjd0 

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
    if self.coords.is_empty() or self.mjd0 is None:
      print('Querying TNS for RA, Dec, and discovery date...')
      json_data = query_tns(tnsname, api_key, tns_id, bot_name)

      if self.coords.is_empty():
        self.coords = Coordinates(json_data['data']['reply']['ra'], json_data['data']['reply']['dec'])
        #print(f'Setting coordinates to TNS coordinates: {self.coords}')
      
      if self.mjd0 is None:
        disc_date = json_data['data']['reply']['discoverydate']
        date = list(disc_date.partition(' '))[0]
        time = list(disc_date.partition(' '))[2]
        date_object = Time(date+"T"+time, format='isot', scale='utc')
        self.mjd0 = date_object.mjd - DISC_DATE_BUFFER
        #print(f'Setting MJD0 to TNS discovery date minus {DISC_DATE_BUFFER}: {self.mjd0}')

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
        result = query_atlas(headers, self.coords.ra.degrees, self.coords.dec.degrees, min_mjd, max_mjd)
        break
      except Exception as e:
        print('Exception caught: '+str(e))
        print('Trying again in 20 seconds! Waiting...')
        time.sleep(20)
        continue
    self.t = result

  # divide the light curve by filter and save into separate files
  def save(self, output_dir, tnsname, overwrite=False):
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

    for filt in ['o','c']:
      filename = get_filename(output_dir, tnsname, filt=filt, control_index=self.control_index)
      indices = lc.ix_equal(colnames=['F'], val=filt)
      print(f'Saving downloaded light curve with filter {filt} (length {len(indices)}) at {filename}...')
      lc.save_by_filename(filename, indices=indices, overwrite=overwrite)

  def __str__(self):
    return f'Full light curve at {self.coords}: control ID = {self.control_index}, MJD0 = {self.mjd0}'

# contains either o-band or c-band measurements
class LightCurve(pdastrostatsclass):
  def __init__(self, control_index=0, filt='o'):
    pdastrostatsclass.__init__(self)
    self.control_index = control_index
    self.filt = filt

  def set_df(self, t:pd.DataFrame):
    self.t = copy.deepcopy(t)

  def check_column_names(self):
    if self.t is None:
      return
    
    for column_name in REQUIRED_COLUMN_NAMES:
      if not column_name in self.t.columns:
        raise RuntimeError(f'ERROR: Missing required column: {column_name}')

  def load(self, output_dir, tnsname):
    filename = get_filename(output_dir, tnsname, self.filt, self.control_index)
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

  def load(self, output_dir, tnsname):
    filename = get_filename(output_dir, tnsname, self.filt, self.control_index, self.mjdbinsize)
    self.load_by_filename(filename)
  
  def save(self, output_dir, tnsname, indices=None, overwrite=False):
    filename = get_filename(output_dir, tnsname, self.filt, self.control_index, self.mjdbinsize)
    self.save_by_filename(filename, indices=indices, overwrite=overwrite)