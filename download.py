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

import requests, argparse, json
import pandas as pd
from getpass import getpass
from lightcurve import Coordinates, SnInfo, FullLightCurve

# define command line arguments
def define_args(parser=None, usage=None, conflict_handler='resolve'):
  if parser is None:
    parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
    
  parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')
  parser.add_argument('--config_file', default='config.ini', type=str, help='file name of ini file with settings for this class')
  parser.add_argument('-l', '--lookbacktime', default=None, type=int, help='lookback time (MJD)')
  parser.add_argument('--min_mjd', default=None, type=float, help='minimum MJD to download')
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
    self.control_coords = None
    self.overwrite = args.overwrite

    config = self.load_config(args)

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
    print(f'Downloading control light curves: {self.controls}')
    if self.controls:
      if args.ctrl_coords:
        pass
      else:
        self.radius = args.radius if args.radius else float(config["download"]["radius"])
        self.num_controls = args.num_controls if args.num_controls else int(config["download"]["num_controls"])
        print(f'Getting circle pattern of {self.num_controls} control light curves with radius of {self.radius}\"')

    """if args.coords:
      parsed_coords = [coord.strip() for coord in args.coords.split(",")]
      if len(parsed_coords) > 2:
        raise RuntimeError('ERROR: Too many coordinates in --coords argument! Please provide comma-separated RA and Dec onlyy.')
      if len(parsed_coords) < 2:
        raise RuntimeError('ERROR: Too few coordinates in --coords argument! Please provide comma-separated RA and Dec.')
      self.settings["coords"] = Coordinates(parsed_coords[0], parsed_coords[1])
    
    if args.mjd0:
      self.settings["mjd0"] = float(args.mjd0)

    if args.lookbacktime:
      self.settings["lookbacktime"] = float(args.lookbacktime)

    # control light curves
    self.controls = False
    if bool(args.controls):
      self.controls = True
      print(f'Downloading control light curves: {self.controls}')
      if args.ctrl_coords:
        self.settings['ctrl_coords'] = args.ctrl_coords
      else:
        self.settings['num_controls'] = int(config["download"]["num_controls"])
        self.settings['radius'] = float(config["download"]["radius"])
        self.settings['flux2mag_sigmalimit'] = float(config["download"]["flux2mag_sigmalimit"])

        # TODO: closebright arguments
    else:
      if args.ctrl_coords:
        print(f'WARNING: In order to download control light curves using the coordinates in {args.ctrl_coords}, use the -c argument')
"""
  def load_config(self, args):
    with open(args.config_file) as config_file:
      return json.load(config_file)
    
  def get_arg_coords(self, args):
    parsed_coords = [coord.strip() for coord in args.coords.split(",")]
    if len(parsed_coords) > 2:
      raise RuntimeError('ERROR: Too many coordinates in --coords argument! Please provide comma-separated RA and Dec onlyy.')
    if len(parsed_coords) < 2:
      raise RuntimeError('ERROR: Too few coordinates in --coords argument! Please provide comma-separated RA and Dec.')
    return Coordinates(parsed_coords[0], parsed_coords[1])

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

  def loop(self, args):
    print('Connecting to ATLAS API...')
    headers = self.connect_atlas()
    if headers is None: 
      raise RuntimeError('ERROR: No token header!')
    
if __name__ == "__main__":
  args = define_args().parse_args()
  download = DownloadLoop(args)
  download.loop(args)