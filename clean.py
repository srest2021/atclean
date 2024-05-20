from typing import Dict, Type, Any
import os, sys, argparse
import pandas as pd
import numpy as np
from lightcurve import SnInfoTable
from download import load_config

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
    output = f'Cut'
    if self.column:
      output += f' on column {self.column}'
    if self.flag:
      output += f' (flag {hex(self.flag)})'
    if self.min_value and self.max_value:
      output += f': min value {self.min_value}, max value {self.max_value}'
    elif self.min_value:
      output += f': min value {self.min_value}'
    elif self.max_value:
      output += f': max value {self.max_value}'
    return output
  
class CutList:
  def __init__(self):
    self.list: Dict[str, Type[Cut]] = {}

  def add_cut(self, cut:Cut, name:str):
    self.list[name] = cut

  def get_cut(self, name:str):
    return self.list[name]

  def can_apply_directly(self, name:str):
    return self.list[name].can_apply_directly()
  
  def __str__(self):
    output = ''
    for cut in self.list:
      output += f'{cut}'

# define command line arguments
def define_args(parser=None, usage=None, conflict_handler='resolve'):
  if parser is None:
    parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
    
  parser.add_argument('tnsnames', nargs='+', help='TNS names of the transients to clean')
  parser.add_argument('--sninfo_file', default=None, type=str, help='file name of .txt file with SN info table')
  parser.add_argument('--config_file', default='config.ini', type=str, help='file name of .ini file with settings for this class')
  parser.add_argument('-o','--overwrite', default=False, action='store_true', help='overwrite existing file with same file name')

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

class CleanLoop:
  def __init__(self, args):
    self.sn = None

    self.tnsnames = args.tnsnames
    print(f'List of transients to clean: {self.tnsnames}')
    if len(self.tnsnames) < 1:
      raise RuntimeError('ERROR: Please specify at least one TNS name to clean.')
    
    self.overwrite = args.overwrite
    print(f'Overwrite existing files: {self.overwrite}')

    config = load_config(args.config_file)

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
    print()
    sninfo_filename = args.sninfo_file if args.sninfo_file else config["dir"]["sninfo_filename"]
    self.sninfo = SnInfoTable(self.input_dir, filename=sninfo_filename)
    
    # TNS credentials
    self.credentials = config["credentials"]
    print(f'TNS ID: {self.credentials["tns_id"]}')
    print(f'TNS bot name: {self.credentials["tns_bot_name"]}')

    # control light curves
    self.controls = bool(args.controls)
    print(f'\nClean control light curves: {self.controls}')
    if self.controls:
      self.num_controls = args.num_controls if args.num_controls else int(config["download"]["num_controls"])
    elif args.num_controls:
      raise RuntimeError('ERROR: Please specify control light curve cleaning (-c or --controls) before using the --num_controls argument.')

    # cuts
    self.cut_list = CutList()
    self.apply_template_correction = args.template_correction
    self.apply_uncert_est = args.uncert_est
    self.apply_uncert_cut = args.uncert_cut
    self.apply_x2_cut = args.x2_cut
    self.apply_controls_cut = args.controls_cut
    self.apply_averaging = args.averaging

    if self.apply_uncert_est:
      pass
      # TODO

    if self.apply_uncert_cut:
      uncert_cut = Cut(column='duJy', 
                       max_value=float(config['uncert_cut']['max_value']), 
                       flag=hexstring_to_int(config['uncert_cut']['flag']))
      self.cut_list.add_cut(uncert_cut, 'uncert_cut')

    if self.apply_x2_cut:
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
      self.cut_list.add_cut(x2_cut, 'x2_cut')
    
    if self.apply_controls_cut:
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
      self.cut_list.add_cut(controls_cut, 'controls_cut')
    
    if self.apply_averaging:
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
      self.cut_list.add_cut(badday_cut, 'badday_cut')

    # TODO
    
  def clean_lcs(self, args, tnsname):
    print(f'\nDOWNLOADING ATLAS LIGHT CURVES FOR: SN {tnsname}')
    # TODO

  def loop(self, args):
    for obj_index in range(len(args.tnsnames)):
      self.clean_lcs(args, args.tnsnames[obj_index]) 

if __name__ == "__main__":
  args = define_args().parse_args()
  clean = CleanLoop(args)
  clean.loop(args)