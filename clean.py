#!/usr/bin/env python

from typing import List, Dict, Type, Any
import os, sys, argparse
import pandas as pd
import numpy as np
from lightcurve import SnInfoTable, Supernova
from download import load_config, make_dir_if_not_exists, parse_comma_separated_string

DEFAULT_CUT_NAMES = ['uncert_cut', 'x2_cut', 'controls_cut', 'badday_cut', 'averaging']

def hexstring_to_int(hexstring):
  return int(hexstring, 16)

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
    if not self.flag or not self.column:
      return False
    if not self.min_value and not self.max_value:
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
  
class CutList:
  def __init__(self):
    self.list: Dict[str, Type[Cut]] = {}

  def add_cut(self, cut:Cut, name:str):
    self.list[name] = cut

  def get_cut(self, name:str):
    return self.list[name]
  
  def has_cut(self, name:str):
    return name in self.list

  def can_apply_directly(self, name:str):
    return self.list[name].can_apply_directly()
  
  def __str__(self):
    output = ''
    for name in self.list:
      output += f'\n{name}: ' + self.list[name].__str__()
    return output

class UncertEstTable():
  def __init__(self, directory, filename=None):
    if filename is None:
      self.filename = f'{directory}/uncert_est_info.txt'
    else:
      self.filename = f'{directory}/{filename}'

    try:
      print(f'\nLoading true uncertainties estimation table at {self.filename}...')
      self.t = pd.read_table(self.filename, delim_whitespace=True)
      print('Success')
    except:
      print(f'No existing true uncertainties estimation table; creating blank table...')
      self.t = pd.DataFrame(columns=['tnsname', 'filter', 'sigma_extra', 'sigma_typical_old', 'sigma_typical_new', 'sigma_typical_new_pct_greater', 'recommended', 'applied'])
    
  def add_row(self, row):
    tnsname = row['tnsname']
    filt = row['filter']

    if len(self.t) > 0:
      matching_ix = np.where(self.t['tnsname'].eq(tnsname) & self.t['filter'].eq(filt))[0]
      if len(matching_ix) > 1:
        raise RuntimeError(f'ERROR: true uncertainties estimation table has {len(matching_ix)} matching rows for TNS name {tnsname} and filter {filt}')
    
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
    print(f'\nSaving true uncertainties estimation table at {self.filename}...')
    self.t.to_string(self.filename)

class ChiSquareCutTable():
  def __init__(self, directory, filename=None):
    if filename is None:
      self.filename = f'{directory}/x2_cut_info.txt'
    else:
      self.filename = f'{directory}/{filename}'

    try:
      print(f'\nLoading chi-square cut table at {self.filename}...')
      self.t = pd.read_table(self.filename, delim_whitespace=True)
      print('Success')
    except:
      print(f'No existing chi-square cut table; creating blank table...')
      self.t = pd.DataFrame(columns=['tnsname', 'filter', 'x2_cut', 'use_preSN_lc', 'stn_bound', 'pct_contamination', 'pct_loss'])

  def add_row(self, row):
    tnsname = row['tnsname']
    filt = row['filter']

    if len(self.t) > 0:
      matching_ix = np.where(self.t['tnsname'].eq(tnsname) & self.t['filter'].eq(filt))[0]
      if len(matching_ix) > 1:
        raise RuntimeError(f'ERROR: chi-square cut table has {len(matching_ix)} matching rows for TNS name {tnsname} and filter {filt}')
    
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
    print(f'\nSaving chi-square cut table at {self.filename}...')
    self.t.to_string(self.filename)

class CleanLoop:
  def __init__(self, 
               input_dir, 
               output_dir, 
               credentials,
               cut_list:CutList, 
               sninfo_filename=None, 
               overwrite=False):
    self.credentials:Dict[str,str] = None
    self.input_dir:str = input_dir
    self.output_dir:str = output_dir
    self.overwrite:bool = overwrite

    self.sn:Supernova = None

    self.cut_list:CutList = cut_list

    self.sninfo:SnInfoTable = SnInfoTable(self.output_dir, filename=sninfo_filename)
    self.uncert_est_info:UncertEstTable = UncertEstTable(self.output_dir)
    if cut_list.has_cut('x2_cut'):
      self.x2_cut_info:ChiSquareCutTable = ChiSquareCutTable(self.output_dir)

  def clean_lcs(self, tnsname, filt, num_controls=0, mjd0=None):
    print(f'\nCLEANING LIGHT CURVES FOR: SN {tnsname}')
    # TODO

    # load the SN light curve, SN info, and control light curves
    self.sn = Supernova(tnsname=tnsname, mjd0=mjd0, filt=filt)
    self.sn.get_tns_data(self.credentials['tns_api_key'], self.credentials['tns_id'], self.credentials['bot_name'])
    self.sn.load_all(self.input_dir, num_controls=num_controls)

  def loop(self, tnsnames, num_controls=0, mjd0=None, filters=['o','c']):
    for obj_index in range(len(tnsnames)):
      for filt in filters:
        self.clean_lcs(tnsnames[obj_index], filt, num_controls=num_controls, mjd0=mjd0) 

def parse_config_filters(args, config):
  if args.filters:
    return parse_comma_separated_string(args.filters)
  else:
    return parse_comma_separated_string(config['convert']['filters'])

def find_config_custom_cuts(config):
  print('\nSearching config file for custom cuts...')
  custom_cuts = []
  for key in config:
    if key.endswith('_cut') and not key in DEFAULT_CUT_NAMES:
      custom_cuts.append(config[key])
  print(f'Found {len(custom_cuts)}')
  return custom_cuts

def parse_config_cuts(args, config):
  cut_list = CutList()
  # TODO
  return cut_list

# define command line arguments
def define_args(parser=None, usage=None, conflict_handler='resolve'):
  if parser is None:
    parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
    
  parser.add_argument('tnsnames', nargs='+', help='TNS names of the transients to clean')
  parser.add_argument('--sninfo_file', default=None, type=str, help='file name of .txt file with SN info table')
  parser.add_argument('--config_file', default='config.ini', type=str, help='file name of .ini file with settings for this class')
  parser.add_argument('-o','--overwrite', default=False, action='store_true', help='overwrite existing file with same file name')
  parser.add_argument('--filters', type=str, default=None, help='comma-separated list of filters to clean')

  # cleaning a single SN and/or controls
  parser.add_argument('--mjd0', type=str, default=None, help='transient start date in MJD')

  # cleaning control light curves
  parser.add_argument('-c','--controls', default=False, action='store_true', help='clean control light curves in addition to transient light curve')
  parser.add_argument('--num_controls', type=int, default=None, help='number of control light curves to load and clean')
  
  # possible cuts
  parser.add_argument('-t', '--template_correction', default=False, action='store_true', help='apply automatic ATLAS template change correction')
  parser.add_argument('-e', '--uncert_est', default=False, action='store_true', help='apply true uncertainty estimation')
  parser.add_argument('-u', '--uncert_cut', default=False, action='store_true', help='apply uncertainty cut')
  parser.add_argument('-x', '--x2_cut', default=False, action='store_true', help='apply chi-square cut')
  parser.add_argument('-n', '--controls_cut', default=False, action='store_true', help='apply control light curve cut')
  parser.add_argument('-g', '--averaging', default=False, action='store_true', help='average light curves and cut bad days')
  parser.add_argument('-m', '--mjd_bin_size', type=float, default=None, help='MJD bin size in days for averaging')
  parser.add_argument('--custom_cuts', default=False, action='store_true', help='scan config file for custom cuts')

  return parser

# TODO: Rewrite this whole thing

"""
# define command line arguments
def define_args(parser=None, usage=None, conflict_handler='resolve'):
  if parser is None:
    parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
    
  parser.add_argument('tnsnames', nargs='+', help='TNS names of the transients to clean')
  parser.add_argument('--sninfo_file', default=None, type=str, help='file name of .txt file with SN info table')
  parser.add_argument('--config_file', default='config.ini', type=str, help='file name of .ini file with settings for this class')
  parser.add_argument('-o','--overwrite', default=False, action='store_true', help='overwrite existing file with same file name')

  # cleaning a single SN and/or controls
  parser.add_argument('--mjd0', type=str, default=None, help='transient start date in MJD')

  # cleaning control light curves
  parser.add_argument('-c','--controls', default=False, action='store_true', help='clean control light curves in addition to transient light curve')
  parser.add_argument('--num_controls', type=int, default=None, help='number of control light curves to load and clean')
  
  # possible cuts
  parser.add_argument('-t', '--template_correction', default=False, action='store_true', help='apply automatic ATLAS template change correction')
  parser.add_argument('-e', '--uncert_est', default=False, action='store_true', help='apply true uncertainty estimation')
  parser.add_argument('-u', '--uncert_cut', default=False, action='store_true', help='apply uncertainty cut')
  parser.add_argument('-x', '--x2_cut', default=False, action='store_true', help='apply chi-square cut')
  parser.add_argument('-c', '--controls_cut', default=False, action='store_true', help='apply control light curve cut')
  parser.add_argument('-g', '--averaging', default=False, action='store_true', help='average light curves and cut bad days')
  parser.add_argument('-m', '--mjd_bin_size', type=float, default=None, help='MJD bin size in days for averaging')
  parser.add_argument('--custom_cuts', default=False, action='store_true', help='scan config file for custom cuts')

  return parser

def parse_settings(args, config):
  settings = {
    'overwrite': args.overwrite,
    'input_dir': config["dir"]["atclean_input"],
    'output_dir': config["dir"]["output"],
    'tns_api_key': config["credentials"]["tns_api_key"],
    'tns_id': config["credentials"]["tns_id"],
    'tns_bot_name': config["credentials"]["tns_bot_name"],
    'sninfo_filename': args.sninfo_file if args.sninfo_file else config["dir"]["sninfo_filename"]
  }

  if len(args.tnsnames) < 1:
    raise RuntimeError('ERROR: Please specify at least one TNS name to clean.')
  settings['tnsnames'] = args.tnsnames
  print(f'\nList of transients to clean: {settings["tnsnames"]}')
  
  print(f'\nOverwrite existing files: {settings["overwrite"]}')
  make_dir_if_not_exists(settings["input_dir"])
  make_dir_if_not_exists(settings["output_dir"])
  print(f'ATClean input directory: {settings["input_dir"]}')
  print(f'Output directory: {settings["output_dir"]}')
  print(f'TNS ID: {settings["tns_id"]}')
  print(f'TNS bot name: {settings["tns_bot_name"]}')

  # control light curves
  settings['apply_to_controls'] = args.controls
  print(f'\nClean control light curves: {settings["apply_to_controls"]}')
  if settings['apply_to_controls']:
    settings['num_controls'] = args.num_controls if args.num_controls else int(config["download"]["num_controls"])
    print(f'Number of control light curves: {settings["num_controls"]}')
  elif args.num_controls:
    raise RuntimeError('ERROR: Please specify control light curve cleaning (-c or --controls) before using the --num_controls argument.')

  settings['apply_template_correction'] = args.template_correction
  print(f'Apply template correction: {settings["apply_template_correction"]}')
  
  settings['apply_uncert_est'] = args.uncert_est
  output = f'Apply true uncertainties estimation: {settings["apply_uncert_est"]}'
  if args.uncert_est: 
    settings['prelim_x2_max_value'] = config['uncert_est']['prelim_x2_max_value']
    output += f' (preliminary chi-square cut: {settings["prelim_x2_max_value"]})'
  print(output)
  
  cut_list = CutList()

  if args.uncert_cut:
    uncert_cut = Cut(column='duJy', 
                     max_value=float(config['uncert_cut']['max_value']), 
                     flag=hexstring_to_int(config['uncert_cut']['flag']))
    cut_list.add_cut(uncert_cut, 'uncert_cut')

  if args.x2_cut:
    params = {
      'stn_cut': float(config['x2_cut']['stn_bound']),
      'cut_start': int(config['x2_cut']['min_cut']),
      'cut_stop': int(config['x2_cut']['max_cut']),
      'cut_step': int(config['x2_cut']['cut_step']),
      'use_preSN_lc': config['x2_cut']['use_preSN_lc'] == 'True'
    }
    x2_cut = Cut(column='chi/N', 
                 max_value=float(config['x2_cut']['max_value']), 
                 flag=hexstring_to_int(config['x2_cut']['flag']), 
                 params=params)
    cut_list.add_cut(x2_cut, 'x2_cut')
  
  if args.controls_cut:
    params = {
      'questionable_flag': hexstring_to_int(config['controls_cut']['questionable_flag']),
      'x2_max': float(config['controls_cut']['x2_max']),
      'x2_flag': hexstring_to_int(config['controls_cut']['x2_flag']),
      'stn_max': float(config['controls_cut']['stn_max']),
      'stn_flag': hexstring_to_int(config['controls_cut']['stn_flag']),
      'Nclip_max': int(config['controls_cut']['Nclip_max']),
      'Nclip_flag': hexstring_to_int(config['controls_cut']['Nclip_flag']),
      'Ngood_min': int(config['controls_cut']['Ngood_min']),
      'Ngood_flag': hexstring_to_int(config['controls_cut']['Ngood_flag'])
    }
    controls_cut = Cut(flag=hexstring_to_int(config['controls_cut']['bad_flag']), 
                       params=params)
    cut_list.add_cut(controls_cut, 'controls_cut')
  
  if args.averaging:
    params = {
      'mjd_bin_size': float(config['averaging']['mjd_bin_size']),
      'x2_max': float(config['averaging']['x2_max']),
      'Nclip_max': int(config['averaging']['Nclip_max']),
      'Ngood_min': int(config['averaging']['Ngood_min']),
      'ixclip_flag': hexstring_to_int(config['averaging']['ixclip_flag']),
      'smallnum_flag': hexstring_to_int(config['averaging']['smallnum_flag'])
    }
    badday_cut = Cut(flag=hexstring_to_int(config['averaging']['flag']),
                     params=params)
    cut_list.add_cut(badday_cut, 'badday_cut')

  if args.custom_cuts:
    custom_cuts_settings = find_custom_cuts_settings(config)
    for i in range(len(custom_cuts_settings)):
      cut_settings = custom_cuts_settings[i]
      try:
        custom_cut = Cut(column=cut_settings['column'], 
                        flag=hexstring_to_int(cut_settings['flag']), 
                        min_value = cut_settings['min_value'] if cut_settings['min_value'] != 'None' else None, 
                        max_value = cut_settings['max_value'] if cut_settings['max_value'] != 'None' else None)
        cut_list.add_cut(custom_cut, f'custom_cut_{i}')
      except Exception as e:
        print(f'WARNING: Could not parse custom cut {cut_settings}: {str(e)}')
  
  print(cut_list,'\n')
  return settings, cut_list

if __name__ == "__main__":
  args = define_args().parse_args()
  config = load_config(args.config_file)
  settings, cut_list = parse_settings(args, config)
  
  clean = CleanLoop(settings, cut_list)
  clean.loop(args.tnsnames)
"""

if __name__ == "__main__":
  args = define_args().parse_args()
  config = load_config(args.config_file)
  
  if len(args.tnsnames) < 1:
    raise RuntimeError('ERROR: Please specify at least one TNS name to clean.')
  print(f'\nList of transients to clean: {args.tnsnames}')

  input_dir = config['dir']['atclean_input']
  output_dir = config['dir']['output']
  sninfo_filename = config['dir']['sninfo_filename']
  make_dir_if_not_exists(input_dir)
  make_dir_if_not_exists(output_dir)
  print(f'\nATClean input directory: {input_dir}')
  print(f'Output directory: {output_dir}')

  print(f'TNS ID: {config["credentials"]["tns_id"]}')
  print(f'TNS bot name: {config["credentials"]["tns_bot_name"]}')

  filters = parse_config_filters(args, config)
  print(f'\nOverwrite existing files: {args.overwrite}')
  print(f'Filters: {filters}')

  # TODO: template correction, uncertainty estimation

  cut_list = parse_config_cuts(args, config)

  print(f'\nClean control light curves: {args.controls}')
  num_controls = 0
  if args.controls:
    num_controls = args.num_controls if args.num_controls else int(config["download"]["num_controls"])
    print(f'Number of control light curves to clean: {num_controls}')
  elif args.num_controls:
    raise RuntimeError('ERROR: Please specify control light curve cleaning (-c or --controls) before using the --num_controls argument.')

  clean = CleanLoop(input_dir, 
                    output_dir, 
                    config['credentials'],
                    cut_list, 
                    sninfo_filename=sninfo_filename, 
                    overwrite=args.overwrite)
  clean.loop(args.tnsnames, 
             num_controls=num_controls,
             mjd0=args.mjd0,
             filters=filters)