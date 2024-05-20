import os, sys, argparse
import pandas as pd
import numpy as np
from lightcurve import LightCurve, AveragedLightCurve

# define command line arguments
def define_args(parser=None, usage=None, conflict_handler='resolve'):
  if parser is None:
    parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
    
  #parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')

  return parser

class CleanLoop:
  def __init__(self, args):
    pass

  def loop(self, args):
    pass

if __name__ == "__main__":
	args = define_args().parse_args()
	clean = CleanLoop(args)
	clean.loop()