#!/usr/bin/env python

import pandas as pd
import numpy as np
from pdastro import pdastrostatsclass, AandB, AnotB
from atlas_lc import atlas_lc

class detect_atlas_lc_bumps():
    def __init__(self, avglc, filt):
        self.avglc = avglc
        self.filt = filt

    	# add simulated bump if necessary and apply rolling gaussian weighted sum to light curve
	def apply_gaussian(self, control_index=0, simparams=None, filt=None):
		bad_ix = avglc.lcs[control_index].ix_unmasked('Mask',self.flags['avg_badday'])

		# make sure there are no lingering simulations
		dropcols=[]
		for col in ['uJysim','SNRsim','simLC','SNRsimsum']:
			if col in avglc.lcs[control_index].t.columns:
				dropcols.append(col)
		if len(dropcols) > 0:
			avglc.lcs[control_index].t.drop(columns=dropcols,inplace=True)

		avglc.lcs[control_index].t['SNR'] = 0.0
		avglc.lcs[control_index].t.loc[bad_ix,'SNR'] = avglc.lcs[control_index].t.loc[bad_ix,'uJy']/avglc.lcs[control_index].t.loc[bad_ix,'duJy']

		if not(simparams is None):
			peakMJDs = simparams['sim_peakMJD'].split(',')
			
			# get the simulated gaussian
			mjds = avglc.lcs[control_index].t.loc[bad_ix,'MJD']
			avglc.lcs[control_index].t.loc[bad_ix,'uJysim'] = avglc.lcs[control_index].t.loc[bad_ix,'uJy']
			avglc.lcs[control_index].t['simLC'] = 0.0
			for peakMJD in peakMJDs:
				peakMJD = float(peakMJD)
				print(f'## Adding simulated gaussian at peak MJD {peakMJD:0.2f} with apparent magnitude {simparams["sim_appmag"]:0.2f}, sigma- of {simparams["sim_sigma_minus"]:0.2f}, and sigma+ of {simparams["sim_sigma_plus"]:0.2f}')

				# get simulated gaussian flux and add to light curve flux
				simflux = gauss2lc(mjds, peakMJD, simparams['sim_sigma_minus'], simparams['sim_sigma_plus'], app_mag=simparams['sim_appmag'])
				avglc.lcs[control_index].t.loc[bad_ix,'uJysim'] += simflux

				# get the simulated lc for all MJDs
				simflux_all = gauss2lc(avglc.lcs[control_index].t['MJDbin'], peakMJD, simparams['sim_sigma_minus'], simparams['sim_sigma_plus'], app_mag=simparams['sim_appmag'])
				avglc.lcs[control_index].t['simLC'] += simflux_all

			# make sure all bad rows have SNRsim = 0.0 so they have no impact on the rolling SNRsum
			avglc.lcs[control_index].t['SNRsim'] = 0.0
			# include simflux in the SNR
			avglc.lcs[control_index].t.loc[bad_ix,'SNRsim'] = avglc.lcs[control_index].t.loc[bad_ix,'uJysim']/avglc.lcs[control_index].t.loc[bad_ix,'duJy']

		gaussian_sigma = round(self.gaussian_sigma/self.mjd_bin_size)
		windowsize = int(6*gaussian_sigma)
		halfwindowsize = int(windowsize*0.5)+1
		print(f'## Sigma (days): {self.gaussian_sigma:0.2f}; MJD bin size (days): {self.mjd_bin_size:0.2f}; sigma (bins): {gaussian_sigma:0.2f}; window size (bins): {windowsize}')

		# calculate the rolling SNR sum
		
		dataindices = np.array(range(len(avglc.lcs[control_index].t)) + np.full(len(avglc.lcs[control_index].t), halfwindowsize))
		
		temp = pd.Series(np.zeros(len(avglc.lcs[control_index].t) + 2*halfwindowsize), name='SNR', dtype=np.float64)
		temp[dataindices] = avglc.lcs[control_index].t['SNR']
		SNRsum = temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=gaussian_sigma)
		avglc.lcs[control_index].t['SNRsum'] = list(SNRsum.loc[dataindices])
		
		norm_temp = pd.Series(np.zeros(len(avglc.lcs[control_index].t) + 2*halfwindowsize), name='norm', dtype=np.float64)
		norm_temp[np.array(range(len(avglc.lcs[control_index].t)) + np.full(len(avglc.lcs[control_index].t), halfwindowsize))] = np.ones(len(avglc.lcs[control_index].t))
		norm_temp_sum = norm_temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=gaussian_sigma)
		
		avglc.lcs[control_index].t['SNRsumnorm'] = list(SNRsum.loc[dataindices] / norm_temp_sum.loc[dataindices] * max(norm_temp_sum.loc[dataindices]))

		# calculate the rolling SNR sum for SNR with simflux
		if not(simparams is None):
			temp = pd.Series(np.zeros(len(avglc.lcs[control_index].t) + 2*halfwindowsize), name='SNRsim', dtype=np.float64)
			temp[dataindices] = avglc.lcs[control_index].t['SNRsim']
			SNRsimsum = temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=gaussian_sigma)
			avglc.lcs[control_index].t['SNRsimsum'] = list(SNRsimsum.loc[dataindices])

		return avglc

    def initialize():
        

if __name__ == "__main__":
	detect_atlas_lc_bumps = detect_atlas_lc_bumps()
	detect_atlas_lc_bumps.initialize()