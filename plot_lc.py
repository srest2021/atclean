#!/usr/bin/env python
"""
@author: Sofia Rest
"""

import pandas as pd
import numpy as np

from light_curve import light_curve
from pdastro import pdastrostatsclass, AnotB

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
#plt.style.use('bmh')

class plot_lc():
	def __init__(self, tnsname, output_dir, flags):
		self.lc = None
		self.filt = None
		self.flags = flags
		self.pdf = PdfPages(f'{output_dir}/{tnsname}/{tnsname}_plots.pdf')

		self.tchange1 = 58417
		self.tchange2 = 58882

	def set(self, lc, filt):
		self.lc = lc 
		self.filt = filt

	def save(self):
		print('\nSaving PDF of plots...')
		self.pdf.close()

	def plot_averaged_lc(self):
		if not self.lc.is_averaged:
			raise RuntimeWarning('ERROR: Light curve to be plotted is not averaged!')
		self.plot_cut_lc(self.flags['flag_badday'])

	def plot_chisquare_cut(self):
		self.plot_cut_lc(self.flags['flag_chisquare'], add2title='chi-square cut')

	def plot_uncertainty_cut(self):
		self.plot_cut_lc(self.flags['flag_uncertainty'], add2title='uncertainty cut')

	def plot_controls_cut(self):
		self.plot_og_control_lcs()
		self.plot_cut_lc(self.flags['flag_controls_bad'], add2title='control light curve cut')

	def plot_all_cuts(self):
		self.plot_cut_lc(self.flags['flag_chisquare']|self.flags['flag_uncertainty']|self.flags['flag_controls_bad'], add2title='all cuts')

	def plot_og_lc(self, add2title=None, xlim_lower=None, xlim_upper=None, ylim_lower=None, ylim_upper=None):
		color = 'orange' if self.filt == 'o' else 'cyan'

		fig = plt.figure(figsize=(10,6), tight_layout=True)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().spines['top'].set_visible(False)
		plt.axhline(linewidth=1,color='k')
		plt.ylabel('Flux (µJy)')
		plt.xlabel('MJD')
		title = f'SN {self.lc.tnsname} {self.filt}-band flux'
		if not(add2title is None):
			title += add2title
		plt.title(title)
		plt.axvline(x=self.tchange1, color='magenta', label='ATLAS template change')
		plt.axvline(x=self.tchange2, color='magenta')

		# set x and y limits
		if xlim_lower is None: xlim_lower = self.lc.pdastro.t['MJD'].min() * 0.999
		if xlim_upper is None: xlim_upper = self.lc.pdastro.t['MJD'].max() * 1.001
		if ylim_lower is None: ylim_lower = self.lc.pdastro.t['uJy'].min()
		if ylim_upper is None: ylim_upper = self.lc.pdastro.t['uJy'].max()
		plt.xlim(xlim_lower,xlim_upper)
		plt.ylim(ylim_lower,ylim_upper)

		plt.errorbar(self.lc.pdastro.t.loc[self.lc.corrected_baseline_ix,'MJD'], self.lc.pdastro.t.loc[self.lc.corrected_baseline_ix,'uJy'], yerr=self.lc.pdastro.t.loc[self.lc.corrected_baseline_ix,'duJy'], fmt='none',ecolor=color,elinewidth=1,c=color)
		plt.scatter(self.lc.pdastro.t.loc[self.lc.corrected_baseline_ix,'MJD'], self.lc.pdastro.t.loc[self.lc.corrected_baseline_ix,'uJy'], s=45,color=color,marker='o',label='Baseline')
		
		plt.errorbar(self.lc.pdastro.t.loc[self.lc.during_sn_ix,'MJD'], self.lc.pdastro.t.loc[self.lc.during_sn_ix,'uJy'], self.lc.pdastro.t.loc[self.lc.during_sn_ix,'duJy'], fmt='none',ecolor='red',elinewidth=1,c='red')
		plt.scatter(self.lc.pdastro.t.loc[self.lc.during_sn_ix,'MJD'], self.lc.pdastro.t.loc[self.lc.during_sn_ix,'uJy'], s=45,color='red',marker='o',label='During SN')
		
		plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3)

		self.pdf.savefig(fig)

	def plot_og_control_lcs(self, add2title=None, xlim_lower=None, xlim_upper=None, ylim_lower=None, ylim_upper=None):
		color = 'orange' if self.filt == 'o' else 'cyan'

		fig = plt.figure(figsize=(10,6), tight_layout=True)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().spines['top'].set_visible(False)
		plt.axhline(linewidth=1,color='k')
		plt.ylabel('Flux (µJy)')
		plt.xlabel('MJD')
		title = f'SN {self.lc.tnsname} and control light curves {self.filt}-band flux'
		if not(add2title is None):
			title += add2title
		plt.title(title)
		plt.axvline(x=self.tchange1, color='magenta', label='ATLAS template change')
		plt.axvline(x=self.tchange2, color='magenta')

		# set x and y limits
		if xlim_lower is None: xlim_lower = self.lc.pdastro.t['MJD'].min() * 0.999
		if xlim_upper is None: xlim_upper = self.lc.pdastro.t['MJD'].max() * 1.001
		if ylim_lower is None: ylim_lower = self.lc.pdastro.t['uJy'].min()
		if ylim_upper is None: ylim_upper = self.lc.pdastro.t['uJy'].max()
		plt.xlim(xlim_lower,xlim_upper)
		plt.ylim(ylim_lower,ylim_upper)
		
		for control_index in range(1, len(self.lc.lcs)+1):
			plt.errorbar(self.lc.lcs[control_index].t['MJD'], self.lc.lcs[control_index].t['uJy'], yerr=self.lc.lcs[control_index].t['duJy'], fmt='none',ecolor='blue',elinewidth=1,c='blue')
			if control_index == 1:
				plt.scatter(self.lc.lcs[control_index].t['MJD'], self.lc.lcs[control_index].t['uJy'], s=45,color='blue',marker='o',label=f'{len(self.lc.lcs)} control light curves')
			else:
				plt.scatter(self.lc.lcs[control_index].t['MJD'], self.lc.lcs[control_index].t['uJy'], s=45,color='blue',marker='o')

		plt.errorbar(self.lc.pdastro.t.loc[self.lc.corrected_baseline_ix,'MJD'], self.lc.pdastro.t.loc[self.lc.corrected_baseline_ix,'uJy'], yerr=self.lc.pdastro.t.loc[self.lc.corrected_baseline_ix,'duJy'], fmt='none',ecolor=color,elinewidth=1,c=color)
		plt.scatter(self.lc.pdastro.t.loc[self.lc.corrected_baseline_ix,'MJD'], self.lc.pdastro.t.loc[self.lc.corrected_baseline_ix,'uJy'], s=45,color=color,marker='o',label='Baseline')
		
		plt.errorbar(self.lc.pdastro.t.loc[self.lc.during_sn_ix,'MJD'], self.lc.pdastro.t.loc[self.lc.during_sn_ix,'uJy'], self.lc.pdastro.t.loc[self.lc.during_sn_ix,'duJy'], fmt='none',ecolor='red',elinewidth=1,c='red')
		plt.scatter(self.lc.pdastro.t.loc[self.lc.during_sn_ix,'MJD'], self.lc.pdastro.t.loc[self.lc.during_sn_ix,'uJy'], s=45,color='red',marker='o',label='During SN')

		plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)

		self.pdf.savefig(fig)

	def plot_cut_lc(self, flags, add2title=None, xlim_lower=None, xlim_upper=None, ylim_lower=None, ylim_upper=None):
		color = 'orange' if self.filt == 'o' else 'cyan'

		good_ix = self.lc.pdastro.ix_unmasked('Mask',maskval=flags)
		bad_ix = AnotB(self.lc.pdastro.getindices(),good_ix)

		fig, (cut, clean) = plt.subplots(1, 2, figsize=(16, 6.5), tight_layout=True)
		title = f'SN {self.lc.tnsname} {self.filt}-band flux'
		if self.lc.is_averaged:
			title += ', averaged'
		if add2title is None:
			title += f', mask value {flags}'
		else:
			title += f', {add2title}'
		plt.suptitle(title, fontsize=19, y=1)
		
		if ylim_lower is None: ylim_lower = -2000
		if ylim_upper is None: 
			ylim_upper = 3*self.lc.get_xth_percentile_flux(95, indices=self.lc.during_sn_ix)
		if xlim_lower is None: xlim_lower = self.lc.discdate - 100
		if xlim_upper is None: xlim_upper = self.lc.discdate + 800
		cut.set_ylim(ylim_lower, ylim_upper)
		cut.set_xlim(xlim_lower,xlim_upper)
		clean.set_ylim(ylim_lower, ylim_upper)
		clean.set_xlim(xlim_lower,xlim_upper)

		#cut.spines.right.set_visible(False)
		#cut.spines.top.set_visible(False)
		cut.errorbar(self.lc.pdastro.t.loc[good_ix,'MJD'], self.lc.pdastro.t.loc[good_ix,'uJy'], yerr=self.lc.pdastro.t.loc[good_ix,'duJy'], fmt='none',ecolor=color,elinewidth=1,c=color)
		cut.scatter(self.lc.pdastro.t.loc[good_ix,'MJD'], self.lc.pdastro.t.loc[good_ix,'uJy'], s=50,color=color,marker='o',label='Kept measurements')
		cut.errorbar(self.lc.pdastro.t.loc[bad_ix,'MJD'], self.lc.pdastro.t.loc[bad_ix,'uJy'], yerr=self.lc.pdastro.t.loc[bad_ix,'duJy'], fmt='none',mfc='white',ecolor=color,elinewidth=1,c=color)
		cut.scatter(self.lc.pdastro.t.loc[bad_ix,'MJD'], self.lc.pdastro.t.loc[bad_ix,'uJy'], s=50,facecolors='white',edgecolors=color,marker='o',label='Cut measurements')
		cut.set_title('All measurements')
		cut.axhline(linewidth=1,color='k')
		cut.set_xlabel('MJD')
		cut.set_ylabel('Flux (uJy)')

		fig.legend(loc='upper center', bbox_to_anchor=(0.5, 0),ncol=2)

		#clean.spines.right.set_visible(False)
		#clean.spines.top.set_visible(False)
		clean.errorbar(self.lc.pdastro.t.loc[good_ix,'MJD'], self.lc.pdastro.t.loc[good_ix,'uJy'], yerr=self.lc.pdastro.t.loc[good_ix,'duJy'], fmt='none',ecolor=color,elinewidth=1,c=color)
		clean.scatter(self.lc.pdastro.t.loc[good_ix,'MJD'], self.lc.pdastro.t.loc[good_ix,'uJy'], s=50,color=color,marker='o',label='Kept measurements')
		clean.set_title('Kept measurements only')
		clean.axhline(linewidth=1,color='k')
		clean.set_xlabel('MJD')
		clean.set_ylabel('Flux (uJy)')
		clean.set_ylim(ylim_lower, ylim_upper)

		self.pdf.savefig(fig, bbox_inches='tight')
