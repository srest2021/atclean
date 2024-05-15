#!/usr/bin/env python
'''
Code adapted from Qinan Wang and Armin Rest by Sofia Rest
'''

"""
Possible inputs:
- .txt table with TNS names and either (RA, Dec, and discovery dates) or (TNS bot credentials)
- list of TNS names in command line and either (.txt table with TNS names, RA, Dec, and discovery dates) or (TNS credentials)

Outputs:
- downloaded light curve files
- if using TNS credentials, new or updated .txt table with TNS names, RA, Dec, and discovery dates
"""

import requests
import pandas as pd
from lightcurve import SnInfo, FullLightCurve

class DownloadLoop:
  def __init__(self, args):
    pass

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
    
    