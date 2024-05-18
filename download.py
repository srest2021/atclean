#!/usr/bin/env python

"""
Code adapted from Qinan Wang and Armin Rest by Sofia Rest

Possible inputs:
- .txt table with TNS names and either (RA, Dec, and MJD0) or (TNS bot credentials)
- list of TNS names in command line and either (.txt table with TNS names, RA, Dec, and MJD0) or (TNS credentials)

Outputs:
- downloaded light curve files
- if using TNS credentials, new or updated .txt table with TNS names, RA, Dec, and MJD0
"""

import os, sys, requests, argparse, configparser
import pandas as pd
from getpass import getpass
from lightcurve import Coordinates, SnInfo, FullLightCurve

class ControlCoordinates:
  def __init__(self):
    self.num_controls = None
    self.radius = None
    self.t = None

  def read(self, filename:str):
    try:
      print(f'Loading control coordinates table at {filename}...')
      self.t = pd.read_table(filename, delim_whitespace=True)
      if not 'ra' in self.t.columns or not 'dec' in self.t.columns:
        raise RuntimeError('ERROR: Control coordinates table must have "ra" and "dec" columns.')
      print('Success')
    except Exception as e:
        raise RuntimeError(f'ERROR: Could not load control coordinates table at {filename}: {str(e)}')

  def construct(self, center_coords:Coordinates, num_controls:int=None, radius:float=None):
    if num_controls:
      self.num_controls = num_controls
    if radius:
      self.radius = radius
    
    self.t = pd.DataFrame(columns=['control_index','ra','dec'])
    # TODO

  def save(self):
    pass
    # TODO

# define command line arguments
def define_args(parser=None, usage=None, conflict_handler='resolve'):
  if parser is None:
    parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
    
  parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')
  parser.add_argument('--config_file', default='config.ini', type=str, help='file name of ini file with settings for this class')
  parser.add_argument('-l', '--lookbacktime', default=None, type=int, help='lookback time (MJD)')
  #parser.add_argument('--min_mjd', default=None, type=float, help='minimum MJD to download')
  parser.add_argument('--max_mjd', default=None, type=float, help='maximum MJD to download')
  parser.add_argument('-o','--overwrite', default=False, action='store_true', help='overwrite existing file with same file name')

  # for downloading single SN only
  parser.add_argument('--coords', type=str, default=None, help='comma-separated RA and Dec of SN light curve to download')
  parser.add_argument('--mjd0', type=str, default=None, help='transient start date in MJD')
  
  # for control light curves
  parser.add_argument('-c','--controls', default=False, action='store_true', help='download control light curves in addition to transient light curve')
  parser.add_argument('-n','--num_controls', default=None, type=float, help='number of control light curves per SN')
  parser.add_argument('-r','--radius', default=None, type=float, help='radius of control light curve circle pattern around SN')
  parser.add_argument('--ctrl_coords', type=str, default=None, help='file name of text file in atclean_input containing table of control light curve coordinates')
  #parser.add_argument('--closebright', type=str, default=None, help='comma-separated RA and Dec coordinates of a nearby bright object interfering with the light curve to become center of control light curve circle')
  
  return parser

class DownloadLoop:
  def __init__(self, args):
    self.lcs = {}
    self.ctrl_coords = ControlCoordinates()
    self.overwrite = args.overwrite

    config = self.load_config(args.config_file)

    self.tnsnames = args.tnsnames
    print(f'List of transients to download from ATLAS: {self.tnsnames}')
    if len(self.tnsnames) < 1:
      raise RuntimeError('ERROR: Please specify at least one TNS name to download.')
    if len(self.tnsnames) > 1 and (not args.coords is None or not args.mjd0 is None or not args.ctrl_coords is None):
      raise RuntimeError(f'ERROR: Cannot specify the same coordinates, MJD0, or control coordinates for multiple SNe in the command line. To run a batch with specific coordinates, use a SN info table at {self.settings["dir"]["atclean_input"]}/{self.settings["dir"]["sninfo_filename"]}.')
    self.flux2mag_sigmalimit = float(config["download"]["flux2mag_sigmalimit"])

    # input/output directories
    self.input_dir = config["dir"]["atclean_input"]
    self.output_dir = config["dir"]["output"]
    print(f'ATClean input directory: {self.input_dir}')
    print(f'Output directory: {self.output_dir}')
    if not os.path.isdir(self.input_dir):
      os.makedirs(self.input_dir)
    if not os.path.isdir(self.output_dir):
      os.makedirs(self.output_dir)
    
    # SN info table
    self.sninfo = SnInfo(self.input_dir, filename=config["dir"]["sninfo_filename"])
    
    # ATLAS and TNS credentials
    self.credentials = config["credentials"]
    print(f'ATLAS username: {self.credentials["atlas_username"]}')
    if self.credentials["atlas_password"]  == 'None':
      self.credentials["atlas_password"] = getpass(prompt='Enter ATLAS password: ')
    print(f'TNS ID: {self.credentials["tns_id"]}')
    print(f'TNS bot name: {self.credentials["tns_bot_name"]}')

    # control light curves
    self.controls = bool(args.controls)
    print(f'Download control light curves: {self.controls}')
    if self.controls:
      if args.ctrl_coords:
        self.ctrl_coords.read(args.ctrl_coords)
      else:
        self.ctrl_coords.radius = args.radius if args.radius else float(config["download"]["radius"])
        self.ctrl_coords.num_controls = args.num_controls if args.num_controls else int(config["download"]["num_controls"])
        print(f'Setting circle pattern of {self.ctrl_coords.num_controls} control light curves with radius of {self.ctrl_coords.radius}\"')

  def load_config(self, config_file):
    cfg = configparser.ConfigParser()
    try:
      print(f'Loading config file at {config_file}...')
      cfg.read(config_file)
    except Exception as e:
      raise RuntimeError(f'ERROR: Could not load config file at {config_file}: {str(e)}')
    return cfg

  def connect_atlas(self):
    baseurl = 'https://fallingstar-data.com/forcedphot'
    resp = requests.post(url=f"{baseurl}/api-token-auth/",data={'username':self.credentials['atlas_username'],'password':self.credentials['atlas_password']})
    if resp.status_code == 200:
      token = resp.json()['token']
      print(f'Token: {token}')
      headers = {'Authorization':f'Token {token}','Accept':'application/json'}
    else:
      raise RuntimeError(f'ERROR in connect_atlas(): {resp.status_code}')
    return headers
  
  def split_arg_coords(self, args):
    parsed_coords = [coord.strip() for coord in args.coords.split(",")]
    if len(parsed_coords) > 2:
      raise RuntimeError('ERROR: Too many coordinates in --coords argument! Please provide comma-separated RA and Dec onlyy.')
    if len(parsed_coords) < 2:
      raise RuntimeError('ERROR: Too few coordinates in --coords argument! Please provide comma-separated RA and Dec.')
    return parsed_coords[0], parsed_coords[1]
  
  def construct_full_lc(self, args, tnsname):
    ra, dec, mjd0 = None, None, None
    
    # first try SN info table
    sninfo_index, sninfo_row = self.sninfo.get_row(tnsname)
    if not sninfo_row is None:
      ra, dec, mjd0 = sninfo_row['ra'], sninfo_row['dec'], sninfo_row['mjd0']
    
    # next try command line args
    if args.coords:
      ra, dec = self.split_arg_coords(args)
    if args.mjd0:
      mjd0 = args.mjd0
    
    try:
      self.lcs[0] = FullLightCurve(0, ra, dec, mjd0)
    except Exception as e:
      print(f'WARNING: Could not construct light curve object with RA {ra}, Dec {dec}, and MJD0 {mjd0}: {str(e)}.')
      self.lcs[0] = FullLightCurve(0)
    
    # try to query TNS for any missing data
    self.lcs[0].get_tns_data(tnsname, 
                             self.credentials['tns_api_key'], 
                             self.credentials['tns_id'], 
                             self.credentials['tns_bot_name'])
  
  def download_lcs(self, args, headers, tnsname):
    print(f'\nDOWNLOADING  ATLAS LIGHT CURVES FOR: SN {tnsname}')

    self.lcs = {}
    
    try:
      self.construct_full_lc(args, tnsname)
    except Exception as e:
      print(f'ERROR: Could not construct light curve object: {str(e)}. Skipping to next SN...')
      return
    
    self.ctrl_coords.construct(self.lcs[0].coords)

    # download SN light curves
    self.lcs[0].download(headers, args.lookbacktime, args.max_mjd)

    # download control light curves
    for i in range(1, len(self.ctrl_coords.t)):
      self.lcs[i] = FullLightCurve(i, 
                                   self.ctrl_coords.t['ra'], 
                                   self.ctrl_coords.t['dec'], 
                                   self.ctrl_coords.t['mjd0'])
      self.lcs[i].download(headers, args.lookbacktime, args.max_mjd)

  def loop(self, args):
    print('\nConnecting to ATLAS API...')
    headers = self.connect_atlas()
    if headers is None: 
      raise RuntimeError('ERROR: No token header!')
    
    for obj_index in range(len(args.tnsnames)):
      self.download_lcs(args, headers, args.tnsnames[obj_index])      
    
if __name__ == "__main__":
  args = define_args().parse_args()
  download = DownloadLoop(args)
  download.loop(args)