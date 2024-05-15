from typing import Dict, Type
import re, json, requests
from astropy import units as u
from astropy.coordinates import Angle
from astropy.time import Time
from collections import OrderedDict

# number of days to subtract from TNS discovery date
# to make sure no SN flux before discovery date
DISC_DATE_BUFFER = 20 

class RA:
  def __init__(self, string:str|None=None):
    self.string:str|None = string
    self.degrees:u.degree|None = self.to_degrees()

  def to_degrees(self) -> u.degree|None:
    if self.string is None:
      return None
    
    s = re.compile('\:')
    if isinstance(self.string,str) and s.search(self.string):
      A = Angle(self.string, u.hour)
    else:
      A = Angle(self.string, u.degree)
    return A.degree
  
class Dec:
  def __init__(self, string:str|None=None):
    self.string:str|None = string
    self.degrees:u.degree|None = self.to_degrees()

  def to_degrees(self) -> u.degree|None:
    if self.string is None:
      return None
    
    A = Angle(self.string, u.degree)
    return A.degree
  
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
    return f'{self.ra.degrees:0.8f} {self.dec.degrees:0.8f}'

class Supernova:
  def __init__(self, tnsname:str=None, ra:str=None, dec:str=None, disc_date:float|None=None):
    self.tnsname = tnsname
    self.num_controls = 0
    self.coords:Coordinates = Coordinates(ra,dec)
    self.disc_date = self.set_disc_date(disc_date)

    self.lcs: Dict[int, Type[LightCurve]] = {}
  
  def get(self, control_index=0):
    try:
      return self.lcs[control_index].t
    except: 
      raise RuntimeError(f"ERROR: Cannot get control light curve {control_index}. Num controls set to {self.num_controls} and {len(self.lcs)} lcs in dictionary.")
  
  def set_disc_date(self, disc_date:float):
    self.disc_date = disc_date - DISC_DATE_BUFFER

  def query_tns(self, api_key, tns_id, bot_name):
    if tns_id == "None" or bot_name == "None":
        raise RuntimeError("Cannot query TNS without TNS ID and bot name. Please specify these parameters in settings.ini.")
    
    try:
      url = 'https://www.wis-tns.org/api/get/object'
      json_file = OrderedDict([("objname",self.tnsname), ("objid",""), ("photometry","1"), ("spectra","1")])
      data = { 'api_key':api_key, 'data':json.dumps(json_file) }
      response = requests.post(url, data=data, headers={'User-Agent': 'tns_marker{"tns_id":"%s","type": "bot", "name":"%s"}' % (tns_id, bot_name)})
      json_data = json.loads(response.text,object_pairs_hook=OrderedDict)
      return json_data
    except Exception as e:
      print(json_data['data']['reply'])
      raise RuntimeError('ERROR in get_tns_json_data(): '+str(e))

  def get_tns_data(self, api_key, tns_id, bot_name):
    if self.coords.is_empty() or self.disc_date is None:
      json_data = self.query_tns(api_key, tns_id, bot_name)

      if self.coords.is_empty():
        self.coords = Coordinates(json_data['data']['reply']['ra'], json_data['data']['reply']['dec'])
      
      if self.disc_date is None:
        disc_date = json_data['data']['reply']['discoverydate']
        date = list(disc_date.partition(' '))[0]
        time = list(disc_date.partition(' '))[2]
        dateobject = Time(date+"T"+time, format='isot', scale='utc')
        self.set_disc_date(dateobject.mjd)

class AveragedSupernova(Supernova):
  def __init__(self, tnsname:str=None, ra:str=None, dec:str=None, disc_date:float|None = None, mjdbinsize:float=1.0):
    Supernova.__init__(self, tnsname, ra, dec, disc_date)
    self.mjdbinsize = mjdbinsize

    self.avg_lcs: Dict[int, Type[AveragedLightCurve]] = {}

  def get_avg(self, control_index:int=0):
    try:
      return self.avg_lcs[control_index].t
    except:
      raise RuntimeError(f"Cannot get averaged control light curve {control_index}. Num controls set to {self.num_controls} and {len(self.avg_lcs)} lcs in dictionary.")
  
class LightCurve:
  def __init__(self, ra:str=None, dec:str=None):
    self.t = None
    self.coords = Coordinates(ra,dec)

class AveragedLightCurve(LightCurve):
  def __init__(self, ra:str=None, dec:str=None, mjdbinsize:float=1.0):
    LightCurve.__init__(self, ra, dec)
    self.mjdbinsize = mjdbinsize