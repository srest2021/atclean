#!/usr/bin/env python

"""
Code adapted from Qinan Wang and Armin Rest by Sofia Rest

Possible inputs:
- .txt table with TNS names and either (RA, Dec, and discovery dates) or (TNS bot credentials)
- list of TNS names in command line and either (.txt table with TNS names, RA, Dec, and discovery dates) or (TNS credentials)

Outputs:
- downloaded light curve files
- if using TNS credentials, new or updated .txt table with TNS names, RA, Dec, and discovery dates
"""

import requests, argparse
import pandas as pd
from lightcurve import SnInfo, FullLightCurve

# define command line arguments
def define_args(parser=None, usage=None, conflict_handler='resolve'):
  if parser is None:
    parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
    
  parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')
  parser.add_argument('--coords', type=str, default=None, help='comma-separated RA and Dec of SN light curve to download')
  parser.add_argument('--discdate', type=str, default=None, help='SN discovery date in MJD')
  
  parser.add_argument('-c','--controls', default=False, action='store_true', help='download control light curves in addition to transient light curve')
  parser.add_argument('--closebright', type=str, default=None, help='comma-separated RA and Dec coordinates of a nearby bright object interfering with the light curve to become center of control light curve circle')
  parser.add_argument('--ctrl_coords', type=str, default=None, help='file name of text file in output_dir containing table of control light curve coordinates')

  parser.add_argument('-f','--cfg_filename', default='settings.ini', type=str, help='file name of ini file with settings for this class')
  parser.add_argument('-l', '--lookbacktime_days', default=None, type=int, help='lookback time (days)')
  parser.add_argument('--mjd_min', default=None, type=float, help='minimum MJD to download')
  parser.add_argument('--mjd_max', default=None, type=float, help='maximum MJD to download')
  parser.add_argument('-o','--overwrite', default=False, action='store_true', help='overwrite existing file with same file name')

  return parser

class DownloadLoop:
  def __init__(self, args):
    self.control_coords = None


  def connect_atlas(self):
    baseurl = 'https://fallingstar-data.com/forcedphot'
    resp = requests.post(url=f"{baseurl}/api-token-auth/",data={'username':self.username,'password':self.password})
    if resp.status_code == 200:
      token = resp.json()['token']
      print(f'Token: {token}')
      headers = {'Authorization':f'Token {token}','Accept':'application/json'}
    else:
      raise RuntimeError(f'ERROR in connect_atlas(): {resp.status_code}')
    return headers

  def loop(self):
    print('Connecting to ATLAS API...')
    headers = self.connect_atlas()
    if headers is None: 
      raise RuntimeError('ERROR: No token header!')
    
if __name__ == "__main__":
  args = define_args().parse_args()
  download = DownloadLoop(args)
  download.loop()