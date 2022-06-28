#!/usr/bin/env python
"""
@author: Sofia Rest
"""

import pandas as pd
import numpy as np

from light_curve import light_curve
from pdastro import pdastrostatsclass

# plotting
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pylab as matlib
import warnings
warnings.simplefilter('error', RuntimeWarning)
warnings.filterwarnings("ignore")
from matplotlib.backends.backend_pdf import PdfPages

# plotting styles
plt.rc('axes', titlesize=18)
plt.rc('axes', labelsize=14)
plt.rc('xtick', labelsize=13)
plt.rc('ytick', labelsize=13)
plt.rc('legend', fontsize=13)
plt.rc('font', size=13)
plt.rcParams['font.size'] = 12
plt.style.use('bmh')

class plot_lc():
	def __init__(self, output_dir, flags):
		self.lc = None
		self.filt = None
		self.flags = flags
		self.pdf = PdfPages(f'{output_dir}/{self.lc.tnsname}/{self.lc.tnsname}_plots.pdf')

		self.tchange1 = 58417
		self.tchange2 = 58882

	def set(self, lc, filt):
		self.lc = lc 
		self.filt = filt

	def save(self):
		self.pdf.close()

	def plot_lc(self, filt, xlim_lower=None, xlim_upper=None, ylim_lower=None, ylim_upper=None):
		color = 'orange' if filt == 'o' else 'cyan'

		fig = plt.figure(figsize=(10,6), tight_layout=True)
	    plt.gca().spines['right'].set_visible(False)
	    plt.gca().spines['top'].set_visible(False)
	    plt.axhline(linewidth=1,color='k')
	    plt.ylabel('Flux (µJy)')
	    plt.xlabel('MJD')
	    title = f'SN {self.lc.tnsname} {filt}-band flux'
	    if not(add2title is None):
	        title += add2title
	    plt.title(title)
	    plt.axvline(x=self.tchange1, color='blue', label='ATLAS template change')
	    plt.axvline(x=self.tchange2, color='blue')

	    # set x and y limits
	    if xlim_lower is None: xlim_lower = self.lc.pdasto.t['MJD'].min() * 0.999
	    if xlim_upper is None: xlim_upper = self.lc.pdasto.t['MJD'].max() * 1.001
	    if ylim_lower is None: ylim_lower = self.lc.pdasto.t['uJy'].min()
	    if ylim_upper is None: ylim_upper = self.lc.pdasto.t['uJy'].max()
	    plt.xlim(xlim_lower,xlim_upper)
	    plt.ylim(ylim_lower,ylim_upper)

	    plt.errorbar(self.lc.pdasto.t.loc[self.lc.corrected_baseline_ix,'MJD'], self.lc.pdasto.t.loc[self.lc.corrected_baseline_ix,'uJy'], yerr=self.lc.pdasto.t.loc[self.lc.corrected_baseline_ix,'duJy'], fmt='none',ecolor=color,elinewidth=1,c=color)
	    plt.scatter(self.lc.pdasto.t.loc[self.lc.corrected_baseline_ix,'MJD'], self.lc.pdasto.t.loc[self.lc.corrected_baseline_ix,'uJy'], s=45,color=color,marker='o',label='Baseline')
	    
	    plt.errorbar(self.lc.pdasto.t.loc[self.lc.during_sn_ix,'MJD'], self.lc.pdasto.t.loc[self.lc.during_sn_ix,'uJy'], self.lc.pdasto.t.loc[self.lc.during_sn_ix,'duJy'], fmt='none',ecolor='red',elinewidth=1,c='red')
	    plt.scatter(self.lc.pdasto.t.loc[self.lc.during_sn_ix,'MJD'], self.lc.pdasto.t.loc[self.lc.during_sn_ix,'uJy'], s=45,color='red',marker='o',label='During SN')
	    
	    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3)

		self.pdf.savefig(fig)

	def plot_cut_lc(self, flags, xlim_lower=None, xlim_upper=None, ylim_lower=None, ylim_upper=None):
		color = 'orange' if self.filt == 'o' else 'cyan'

		good_ix = self.lc.pdastro.ix_unmasked('Mask',maskval=flags)
	    bad_ix = AnotB(self.lc.pdastro.getindices(),good_ix)

	    fig, (cut, clean) = plt.subplots(1, 2, figsize=(16, 6.5), tight_layout=True)
	    title = f'SN {self.lc.tnsname} {self.filt}-band'
	    if self.lc.is_averaged:
	        title += ', averaged'
	    title += f', mask value {mask}'
	    plt.suptitle(title, fontsize=19, y=1)
	    
	    if ylim_lower is None: ylim_lower = -2000
	    if ylim_upper is None: 
	        ylim_upper = 3*self.lc.get_xth_percentile_flux(95, indices=self.lc.afterdiscdate_ix)
	    if xlim_lower is None: xlim_lower = self.lc.discdate - 100
	    if xlim_upper is None: xlim_upper = self.lc.discdate + 800
	    cut.set_ylim(ylim_lower, ylim_upper)
	    cut.set_xlim(xlim_lower,xlim_upper)
	    clean.set_ylim(ylim_lower, ylim_upper)
	    clean.set_xlim(xlim_lower,xlim_upper)

	    cut.spines.right.set_visible(False)
	    cut.spines.top.set_visible(False)
	    cut.errorbar(lc.pdastro.t.loc[good_ix,'MJD'], lc.pdastro.t.loc[good_ix,'uJy'], yerr=lc.pdastro.t.loc[good_ix,'duJy'], fmt='none',ecolor=color,elinewidth=1,c=color)
	    cut.scatter(lc.pdastro.t.loc[good_ix,'MJD'], lc.pdastro.t.loc[good_ix,'uJy'], s=50,color=color,marker='o',label='Kept measurements')
	    cut.errorbar(lc.pdastro.t.loc[bad_ix,'MJD'], lc.pdastro.t.loc[bad_ix,'uJy'], yerr=lc.pdastro.t.loc[bad_ix,'duJy'], fmt='none',mfc='white',ecolor=color,elinewidth=1,c=color)
	    cut.scatter(lc.pdastro.t.loc[bad_ix,'MJD'], lc.pdastro.t.loc[bad_ix,'uJy'], s=50,facecolors='white',edgecolors=color,marker='o',label='Cut measurements')
	    cut.set_title('All measurements')
	    cut.axhline(linewidth=1,color='k')
	    cut.set_xlabel('MJD')
	    cut.set_ylabel('Flux (uJy)')

	    fig.legend(loc='upper center', bbox_to_anchor=(0.5, 0),ncol=2)

	    clean.spines.right.set_visible(False)
	    clean.spines.top.set_visible(False)
	    clean.errorbar(lc.pdastro.t.loc[good_ix,'MJD'], lc.pdastro.t.loc[good_ix,'uJy'], yerr=lc.pdastro.t.loc[good_ix,'duJy'], fmt='none',ecolor=color,elinewidth=1,c=color)
	    clean.scatter(lc.pdastro.t.loc[good_ix,'MJD'], lc.pdastro.t.loc[good_ix,'uJy'], s=50,color=color,marker='o',label='Kept measurements')
	    clean.set_title('Kept measurements only')
	    clean.axhline(linewidth=1,color='k')
	    clean.set_xlabel('MJD')
	    clean.set_ylabel('Flux (uJy)')
	    clean.set_ylim(ylim_lower, ylim_upper)

		self.pdf.savefig(fig)

	def plot_averaged_lc(self):
		self.plot_cut_lc(self.flags['flag_badday'])

	def plot_chisquare_cut(self):
		self.plot_cut_lc(self.flags['flag_chisquare'])

	def plot_uncertainty_cut(self):
		self.plot_cut_lc(self.flags['flag_uncertainty'])

	def plot_controls_cut(self):
		self.plot_cut_lc(self.flags['flag_controls_bad'])
