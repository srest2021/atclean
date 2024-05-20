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

from typing import Dict, Type
import os, sys, requests, argparse, configparser, math
import pandas as pd
import numpy as np
from getpass import getpass
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from lightcurve import Coordinates, SnInfoTable, FullLightCurve

"""
UTILITY
"""

def make_dir_if_not_exists(directory):
  if not os.path.isdir(directory):
    os.makedirs(directory)

def load_config(config_file):
    cfg = configparser.ConfigParser()
    try:
      print(f'\nLoading config file at {config_file}...')
      cfg.read(config_file)
    except Exception as e:
      raise RuntimeError(f'ERROR: Could not load config file at {config_file}: {str(e)}')
    return cfg

class ControlCoordinatesTable:
  def __init__(self):
    self.num_controls = None
    self.radius = None
    self.t = None

    self.closebright_min_dist = None

  def read(self, filename:str):
    try:
      print(f'Loading control coordinates table at {filename}...')
      self.t = pd.read_table(filename, delim_whitespace=True)
      if not 'ra' in self.t.columns or not 'dec' in self.t.columns:
        raise RuntimeError('ERROR: Control coordinates table must have "ra" and "dec" columns.')
      print('Success')
    except Exception as e:
        raise RuntimeError(f'ERROR: Could not load control coordinates table at {filename}: {str(e)}')

    self.num_controls = len(self.t)
    self.t['control_index'] = range(1,self.num_controls+1)
    
    with pd.option_context('display.float_format', '{:,.8f}'.format):
      print('Control light curve coordinates read from file: \n',self.t.to_string())

  def update_row(self, control_index:int, full_control_lc:FullLightCurve):
    ix = np.where(self.t['control_index'] == control_index)[0]
    if len(ix) > 1:
      raise RuntimeError(f'ERROR: Cannot update row in control coordinates table for control index {control_index}: duplicate rows.')
    index = ix[0]

    # update corresponding row in table with n_detec, n_detec_o, and n_detec_c counts
    total_len, o_len, c_len = full_control_lc.get_filt_lens()
    self.t[index, 'n_detec'] = total_len
    self.t[index, 'n_detec_o'] = o_len
    self.t[index, 'n_detec_c'] = c_len

  def add_row(self, tnsname:str, control_index:int, coords:Coordinates, ra_offset=0, dec_offset=0, radius=0, n_detec=0, n_detec_o=0, n_detec_c=0):
    row = {
      'tnsname': tnsname,
      'control_index': control_index,
      'ra': f'{coords.ra.angle.degree:0.14f}',
      'dec': f'{coords.dec.angle.degree:0.14f}',
      'ra_offset': f'{ra_offset.degree:0.14f}' if isinstance(ra_offset, Angle) else ra_offset,
      'dec_offset': f'{dec_offset.degree:0.14f}' if isinstance(dec_offset, Angle) else dec_offset,
      'radius_arcsec': radius.arcsecond if isinstance(radius, Angle) else radius,
      'n_detec': n_detec,
      'n_detec_o': n_detec_o,
      'n_detec_c': n_detec_c
    }
    
    self.t = pd.concat([self.t, pd.DataFrame([row])], ignore_index=True)

  def get_distance(self, coord1:Coordinates, coord2:Coordinates) -> Angle:
    c1 = SkyCoord(coord1.ra.angle, coord1.dec.angle, frame='fk5')
    c2 = SkyCoord(coord2.ra.angle, coord2.dec.angle, frame='fk5')
    return c1.separation(c2)

  def construct_row(self, i:int, sn_coords:Coordinates, center_coords:Coordinates, r:Angle, closebright=False):
    angle = Angle(i*360.0 / self.num_controls, u.degree)
    
    ra_distance = Angle(r.degree * math.cos(angle.radian), u.degree)
    ra_offset = Angle(ra_distance.degree * (1.0/math.cos(center_coords.dec.angle.radian)), u.degree)
    ra = Angle(center_coords.ra.angle.degree + ra_offset.degree, u.degree)

    dec_offset = Angle(r.degree * math.sin(angle.radian), u.degree)
    dec = Angle(center_coords.dec.angle.degree + dec_offset.degree, u.degree)

    coords = Coordinates()
    coords.ra.angle = ra
    coords.dec.angle = dec

    if closebright: # check to see if control light curve location is within minimum distance from SN location
      offset_sep = self.get_distance(sn_coords, coords).arcsecond
      if offset_sep < self.closebright_min_dist:
        print(f'Control light curve {i:3d} too close to SN location ({offset_sep}\" away) with minimum distance to SN as {self.closebright_min_dist}\"; skipping control light curve...')
        return
    
    self.add_row(np.nan, i, coords, ra_offset=ra_offset, dec_offset=dec_offset, radius=r)

  def construct(self, full_sn_lc:FullLightCurve, tnsname:str, center_coords:Coordinates, num_controls:int=None, radius:float=None, closebright=False):
    if num_controls:
      self.num_controls = num_controls
    if radius:
      self.radius = radius
    
    self.t = pd.DataFrame(columns=['tnsname','control_index','ra','dec','ra_offset','dec_offset','radius_arcsec','n_detec','n_detec_o','n_detec_c'])

    # add row for SN position
    total_len, o_len, c_len = full_sn_lc.get_filt_lens()
    if closebright:
      # circle pattern radius is distance between SN and bright object
      r = self.get_distance(full_sn_lc.coords, center_coords)

      self.add_row(tnsname, 0, full_sn_lc.coords, ra_offset=np.nan, dec_offset=np.nan, radius=r, n_detec=total_len, n_detec_o=o_len, n_detec_c=c_len)
    else:
      r = Angle(self.radius, u.arcsec)
      self.add_row(tnsname, 0, full_sn_lc.coords, n_detec=total_len, n_detec_o=o_len, n_detec_c=c_len)

    # add row for each control light curve
    for i in range(1,self.num_controls+1):
      self.construct_row(i, full_sn_lc.coords, center_coords, r, closebright=closebright)

    with pd.option_context('display.float_format', '{:,.8f}'.format):
      print('Control light curve coordinates generated: \n',self.t.to_string())

  def save(self, directory, filename=None, tnsname=None):
    if filename is None:
      if tnsname is None:
        raise RuntimeError('ERROR: Please provide either a filename or a TNS name to save the control coordinates table.')
      filename = f'{directory}/{tnsname}/{tnsname}_control_coords.txt'
    else: 
      filename = f'{directory}/{tnsname}/{filename}'

    print(f'Saving control coordinates table at {filename}...')
    self.t.to_string(filename, index=False)

"""
DOWNLOADING ATLAS LIGHT CURVES
"""

# define command line arguments
def define_args(parser=None, usage=None, conflict_handler='resolve'):
  if parser is None:
    parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
    
  parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')
  parser.add_argument('--sninfo_file', default=None, type=str, help='file name of .txt file with SN info table')
  parser.add_argument('--config_file', default='config.ini', type=str, help='file name of .ini file with settings for this class')
  parser.add_argument('-l', '--lookbacktime', default=None, type=int, help='lookback time (MJD)')
  #parser.add_argument('--min_mjd', default=None, type=float, help='minimum MJD to download')
  parser.add_argument('--max_mjd', default=None, type=float, help='maximum MJD to download')
  parser.add_argument('-o','--overwrite', default=False, action='store_true', help='overwrite existing file with same file name')

  # for downloading single SN only
  parser.add_argument('--coords', type=str, default=None, help='comma-separated RA and Dec of SN light curve to download')
  parser.add_argument('--mjd0', type=str, default=None, help='transient start date in MJD')
  
  # for control light curves
  parser.add_argument('-c','--controls', default=False, action='store_true', help='download control light curves in addition to transient light curve')
  parser.add_argument('-n','--num_controls', default=None, type=int, help='number of control light curves per SN')
  parser.add_argument('-r','--radius', default=None, type=float, help='radius of control light curve circle pattern around SN')
  parser.add_argument('--ctrl_coords', type=str, default=None, help='file name of text file in atclean_input containing table of control light curve coordinates')
  # for downloading single SN with control light curves only
  parser.add_argument('--closebright', type=str, default=None, help='comma-separated RA and Dec coordinates of a nearby bright object interfering with the light curve to become center of control light curve circle')
  
  return parser

class DownloadLoop:
  def __init__(self, args):
    self.lcs: Dict[int, Type[FullLightCurve]] = {}
    self.ctrl_coords = ControlCoordinatesTable()

    self.tnsnames = args.tnsnames
    print(f'List of transients to download from ATLAS: {self.tnsnames}')
    if len(self.tnsnames) < 1:
      raise RuntimeError('ERROR: Please specify at least one TNS name to download.')
    if len(self.tnsnames) > 1 and (not args.coords is None or not args.mjd0 is None or not args.ctrl_coords is None or not args.closebright is None):
      raise RuntimeError(f'ERROR: Cannot specify the same coordinates, MJD0, or control/closebright coordinates for multiple SNe in the command line. To run a batch with specific coordinates, use a SN info table at {self.settings["dir"]["atclean_input"]}/{self.settings["dir"]["sninfo_filename"]}.')

    self.overwrite = args.overwrite
    print(f'Overwrite existing files: {self.overwrite}')
    if args.lookbacktime:
      print(f'Lookback time (days): {args.lookbacktime}')
    if args.max_mjd:
      print(f'Max MJD: {args.max_mjd} MJD')

    config = load_config(args.config_file)
    self.flux2mag_sigmalimit = float(config["download"]["flux2mag_sigmalimit"])
    print(f'Sigma limit when converting flux to magnitude : {self.flux2mag_sigmalimit}')
    
    # input/output directories
    self.input_dir = config["dir"]["atclean_input"]
    self.output_dir = config["dir"]["output"]
    print(f'ATClean input directory: {self.input_dir}')
    print(f'Output directory: {self.output_dir}')
    make_dir_if_not_exists(self.input_dir)
    make_dir_if_not_exists(self.output_dir)
    
    # SN info table
    print()
    sninfo_filename = args.sninfo_file if args.sninfo_file else config["dir"]["sninfo_filename"]
    self.sninfo = SnInfoTable(self.output_dir, filename=sninfo_filename)
    
    # ATLAS and TNS credentials
    self.credentials = config["credentials"]
    print(f'\nATLAS username: {self.credentials["atlas_username"]}')
    if self.credentials["atlas_password"]  == 'None':
      self.credentials["atlas_password"] = getpass(prompt='Enter ATLAS password: ')
    print(f'TNS ID: {self.credentials["tns_id"]}')
    print(f'TNS bot name: {self.credentials["tns_bot_name"]}')

    # control light curves
    self.controls = args.controls
    print(f'\nDownload control light curves: {self.controls}')
    if self.controls:
      if args.ctrl_coords:
        self.ctrl_coords.read(args.ctrl_coords)
      else:
        if args.closebright:
          self.ctrl_coords.closebright_min_dist = float(config["download"]["closebright_min_dist"])
          print(f'Circle pattern centered around close bright object at {args.closebright}; minimum distance from SN set to {self.ctrl_coords.closebright_min_dist}')
        self.ctrl_coords.radius = args.radius if args.radius else float(config["download"]["radius"])
        self.ctrl_coords.num_controls = args.num_controls if args.num_controls else int(config["download"]["num_controls"])
        print(f'Setting circle pattern of {self.ctrl_coords.num_controls:d} control light curves', end='')
        if args.closebright:
          print(f'with radius of {self.ctrl_coords.radius}\" from center')
    elif args.ctrl_coords or args.closebright or args.num_controls or args.radius:
      raise RuntimeError('ERROR: Please specify control light curve downloading (-c or --controls) before using any of the following arguments: --ctrl_coords, --closebright, --num_controls, --radius.')

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
  
  def split_arg_coords(self, arg_coords):
    parsed_coords = [coord.strip() for coord in arg_coords.split(",")]
    if len(parsed_coords) > 2:
      raise RuntimeError('ERROR: Too many coordinates in --coords argument! Please provide comma-separated RA and Dec onlyy.')
    if len(parsed_coords) < 2:
      raise RuntimeError('ERROR: Too few coordinates in --coords argument! Please provide comma-separated RA and Dec.')
    return parsed_coords[0], parsed_coords[1]
  
  def construct_full_lc(self, args, tnsname):
    ra, dec, mjd0 = None, None, None
    
    # first try SN info table
    _, sninfo_row = self.sninfo.get_row(tnsname)
    if not sninfo_row is None:
      ra, dec, mjd0 = sninfo_row['ra'], sninfo_row['dec'], sninfo_row['mjd0']

    # next try command line args
    if args.coords:
      ra, dec = self.split_arg_coords(args.coords)
      print(f'Setting coordinates to --coords argument: RA {ra}, Dec {dec}')
    if args.mjd0:
      mjd0 = args.mjd0
      print(f'Setting MJD0 to --mjd0 argument: {mjd0} MJD')
    
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
    
    # add final RA, Dec, MJD0 to SN info table
    self.sninfo.add_row_info(tnsname, self.lcs[0].coords, self.lcs[0].mjd0)
  
  def download_lcs(self, args, headers, tnsname):
    print(f'\nDOWNLOADING ATLAS LIGHT CURVES FOR: SN {tnsname}')

    self.lcs = {}
    try:
      self.construct_full_lc(args, tnsname)
    except Exception as e:
      print(f'ERROR: Could not construct light curve object: {str(e)}. Skipping to next SN...')
      return

    # download SN light curves
    self.lcs[0].download(headers, lookbacktime=args.lookbacktime, max_mjd=args.max_mjd)
    self.lcs[0].save(self.output_dir, tnsname, overwrite=args.overwrite)

    # save SN info table
    self.sninfo.save()

    if args.controls and not args.ctrl_coords:
      # construct control coordinates table
      if args.closebright:
        # TODO: option to parse from SN info table
        parsed_ra, parsed_dec = self.split_arg_coords(args.closebright)
        center_coords = Coordinates(parsed_ra, parsed_dec)
        self.ctrl_coords.construct(self.lcs[0], tnsname, center_coords, closebright=True)
      else:
        self.ctrl_coords.construct(self.lcs[0], tnsname, self.lcs[0].coords, closebright=False)

    if args.controls:
      # download control light curves
      for i in range(1, len(self.ctrl_coords.t)):
        print(f'\nControl light curve {i}')
        self.lcs[i] = FullLightCurve(i, 
                                    self.ctrl_coords.t['ra'], 
                                    self.ctrl_coords.t['dec'], 
                                    self.ctrl_coords.t['mjd0'])
        self.lcs[i].download(headers, lookbacktime=args.lookbacktime, max_mjd=args.max_mjd)
        self.lcs[i].save(self.output_dir, tnsname, overwrite=args.overwrite)
        self.ctrl_coords.update_row(i, self.lcs[i])

      # save control coordinates table
      self.ctrl_coords.save(self.output_dir, tnsname=tnsname)

  def loop(self, args):
    print('\nConnecting to ATLAS API...')
    headers = self.connect_atlas()
    if headers is None: 
      raise RuntimeError('ERROR: No token header!')
    #headers = {}
    
    for obj_index in range(len(args.tnsnames)):
      self.download_lcs(args, headers, args.tnsnames[obj_index])      
    
if __name__ == "__main__":
  args = define_args().parse_args()
  download = DownloadLoop(args)
  download.loop(args)