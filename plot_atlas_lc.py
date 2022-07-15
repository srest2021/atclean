#!/usr/bin/env python
"""
@author: Sofia Rest
"""

import sys
import pandas as pd
import numpy as np

from atlas_lc import atlas_lc
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

class plot_atlas_lc():
	def __init__(self, tnsname, output_dir, args, add2filename=None, flags=None):
		self.lc = None
		self.filt = None
		self.flags = flags
		self.args = args

		if add2filename is None:
			self.pdf = PdfPages(f'{output_dir}/{tnsname}/{tnsname}_plots.pdf')
		else:
			self.pdf = PdfPages(f'{output_dir}/{tnsname}/{tnsname}_plots_{add2filename}.pdf')

		self.tchange1 = 58417
		self.tchange2 = 58882

		self.xlim_lower = None
		self.xlim_upper = None
		self.ylim_lower = None
		self.ylim_upper = None

	def set(self, lc, filt):
		self.lc = lc 
		self.filt = filt

		if not self.lc.is_averaged:
			self.set_plot_lims()

	def save(self):
		print('\nSaving PDF of plots...')
		self.pdf.close()

	def set_plot_lims(self):
		self.xlim_lower = self.lc.discdate-200
		self.xlim_upper = self.lc.discdate+800
		self.ylim_lower = 3*self.lc.get_xth_percentile_flux(1, indices=self.lc.during_sn_ix)
		self.ylim_upper = 3*self.lc.get_xth_percentile_flux(97, indices=self.lc.during_sn_ix)

	def plot_averaged_lc(self):
		if not self.lc.is_averaged:
			raise RuntimeWarning('ERROR: Light curve to be plotted is not averaged!')
		self.plot_cut_lc(self.flags['avg_badday'], add2title=f'MJD bin size {self.lc.mjd_bin_size:0.1f} days')

	def plot_chisquare_cut(self):
		self.plot_cut_lc(self.flags['chisquare'], add2title='chi-square cut')

	def plot_uncertainty_cut(self):
		self.plot_cut_lc(self.flags['uncertainty'], add2title='uncertainty cut')

	def plot_controls_cut(self, num_controls):
		self.plot_og_control_lcs(num_controls)
		self.plot_cut_lc(self.flags['controls_bad'], add2title='control light curve cut')

	def plot_all_cuts(self):
		self.plot_cut_lc(self.flags['chisquare']|self.flags['uncertainty']|self.flags['controls_bad'], add2title='all cuts')

	def plot_og_lc(self, separate_baseline=True, add2title=None): #, xlim_lower=None, xlim_upper=None, ylim_lower=None, ylim_upper=None):
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

		color = 'orange' if self.filt == 'o' else 'cyan'

		if separate_baseline:
			plt.errorbar(self.lc.lcs[0].t.loc[self.lc.corrected_baseline_ix,'MJD'], self.lc.lcs[0].t.loc[self.lc.corrected_baseline_ix,'uJy'], yerr=self.lc.lcs[0].t.loc[self.lc.corrected_baseline_ix,'duJy'], fmt='none',ecolor=color,elinewidth=1,c=color)
			plt.scatter(self.lc.lcs[0].t.loc[self.lc.corrected_baseline_ix,'MJD'], self.lc.lcs[0].t.loc[self.lc.corrected_baseline_ix,'uJy'], s=45,color=color,marker='o',label='Baseline')
			
			plt.errorbar(self.lc.lcs[0].t.loc[self.lc.during_sn_ix,'MJD'], self.lc.lcs[0].t.loc[self.lc.during_sn_ix,'uJy'], self.lc.lcs[0].t.loc[self.lc.during_sn_ix,'duJy'], fmt='none',ecolor='red',elinewidth=1,c='red')
			plt.scatter(self.lc.lcs[0].t.loc[self.lc.during_sn_ix,'MJD'], self.lc.lcs[0].t.loc[self.lc.during_sn_ix,'uJy'], s=45,color='red',marker='o',label='During SN')

			plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3)
		else:
			plt.errorbar(self.lc.lcs[0].t['MJD'], self.lc.lcs[0].t['uJy'], self.lc.lcs[0].t['duJy'], fmt='none',ecolor=color,elinewidth=1,c=color)
			plt.scatter(self.lc.lcs[0].t['MJD'], self.lc.lcs[0].t['uJy'], s=45,color=color,marker='o')

		xlim_lower = self.args.xlim_lower if not(self.args.xlim_lower is None) else self.xlim_lower
		xlim_upper = self.args.xlim_upper if not(self.args.xlim_upper is None) else self.xlim_upper
		plt.xlim(xlim_lower, xlim_upper)
		ylim_lower = self.args.ylim_lower if not(self.args.ylim_lower is None) else self.ylim_lower
		ylim_upper = self.args.ylim_upper if not(self.args.ylim_upper is None) else self.ylim_upper
		plt.ylim(ylim_lower, ylim_upper)

		self.pdf.savefig(fig)

	def plot_og_control_lcs(self, num_controls, add2title=None): #, xlim_lower=None, xlim_upper=None, ylim_lower=None, ylim_upper=None):
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
		
		for control_index in range(1, num_controls+1):
			plt.errorbar(self.lc.lcs[control_index].t['MJD'], self.lc.lcs[control_index].t['uJy'], yerr=self.lc.lcs[control_index].t['duJy'], fmt='none',ecolor='blue',elinewidth=1,c='blue')
			if control_index == 1:
				plt.scatter(self.lc.lcs[control_index].t['MJD'], self.lc.lcs[control_index].t['uJy'], s=45,color='blue',marker='o',label=f'{len(self.lc.lcs)} control light curves')
			else:
				plt.scatter(self.lc.lcs[control_index].t['MJD'], self.lc.lcs[control_index].t['uJy'], s=45,color='blue',marker='o')

		color = 'orange' if self.filt == 'o' else 'cyan'

		plt.errorbar(self.lc.lcs[0].t.loc[self.lc.corrected_baseline_ix,'MJD'], self.lc.lcs[0].t.loc[self.lc.corrected_baseline_ix,'uJy'], yerr=self.lc.lcs[0].t.loc[self.lc.corrected_baseline_ix,'duJy'], fmt='none',ecolor=color,elinewidth=1,c=color)
		plt.scatter(self.lc.lcs[0].t.loc[self.lc.corrected_baseline_ix,'MJD'], self.lc.lcs[0].t.loc[self.lc.corrected_baseline_ix,'uJy'], s=45,color=color,marker='o',label='Baseline')
		
		plt.errorbar(self.lc.lcs[0].t.loc[self.lc.during_sn_ix,'MJD'], self.lc.lcs[0].t.loc[self.lc.during_sn_ix,'uJy'], self.lc.lcs[0].t.loc[self.lc.during_sn_ix,'duJy'], fmt='none',ecolor='red',elinewidth=1,c='red')
		plt.scatter(self.lc.lcs[0].t.loc[self.lc.during_sn_ix,'MJD'], self.lc.lcs[0].t.loc[self.lc.during_sn_ix,'uJy'], s=45,color='red',marker='o',label='During SN')

		plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)

		xlim_lower = self.args.xlim_lower if not(self.args.xlim_lower is None) else self.xlim_lower
		xlim_upper = self.args.xlim_upper if not(self.args.xlim_upper is None) else self.xlim_upper
		plt.xlim(xlim_lower, xlim_upper)
		ylim_lower = self.args.ylim_lower if not(self.args.ylim_lower is None) else self.ylim_lower
		ylim_upper = self.args.ylim_upper if not(self.args.ylim_upper is None) else self.ylim_upper
		plt.ylim(ylim_lower, ylim_upper)

		self.pdf.savefig(fig)

	def plot_cut_lc(self, flags, add2title=None): #, xlim_lower=None, xlim_upper=None, ylim_lower=None, ylim_upper=None):
		good_ix = self.lc.lcs[0].ix_unmasked('Mask',maskval=flags)
		bad_ix = AnotB(self.lc.lcs[0].getindices(),good_ix)

		fig, (cut, clean) = plt.subplots(1, 2, figsize=(16, 6.5), tight_layout=True)
		title = f'SN {self.lc.tnsname} {self.filt}-band flux'
		if self.lc.is_averaged:
			title += ', averaged'
		if add2title is None:
			title += f', mask value {flags}'
		else:
			title += f', {add2title}'
		plt.suptitle(title, fontsize=19, y=1)

		# set x and y limits
		xlim_lower = self.args.xlim_lower if not(self.args.xlim_lower is None) else self.xlim_lower
		xlim_upper = self.args.xlim_upper if not(self.args.xlim_upper is None) else self.xlim_upper
		cut.set_xlim(xlim_lower, xlim_upper)
		clean.set_xlim(xlim_lower, xlim_upper)
		ylim_lower = self.args.ylim_lower if not(self.args.ylim_lower is None) else self.ylim_lower
		ylim_upper = self.args.ylim_upper if not(self.args.ylim_upper is None) else self.ylim_upper
		cut.set_ylim(ylim_lower, ylim_upper)
		clean.set_ylim(ylim_lower, ylim_upper)

		color = 'orange' if self.filt == 'o' else 'cyan'

		cut.errorbar(self.lc.lcs[0].t.loc[good_ix,'MJD'], self.lc.lcs[0].t.loc[good_ix,'uJy'], yerr=self.lc.lcs[0].t.loc[good_ix,'duJy'], fmt='none',ecolor=color,elinewidth=1,c=color)
		cut.scatter(self.lc.lcs[0].t.loc[good_ix,'MJD'], self.lc.lcs[0].t.loc[good_ix,'uJy'], s=50,color=color,marker='o',label='Kept measurements')
		cut.errorbar(self.lc.lcs[0].t.loc[bad_ix,'MJD'], self.lc.lcs[0].t.loc[bad_ix,'uJy'], yerr=self.lc.lcs[0].t.loc[bad_ix,'duJy'], fmt='none',mfc='white',ecolor=color,elinewidth=1,c=color)
		cut.scatter(self.lc.lcs[0].t.loc[bad_ix,'MJD'], self.lc.lcs[0].t.loc[bad_ix,'uJy'], s=50,facecolors='white',edgecolors=color,marker='o',label='Cut measurements')
		cut.set_title('All measurements')
		cut.axhline(linewidth=1,color='k')
		cut.set_xlabel('MJD')
		cut.set_ylabel('Flux (uJy)')

		fig.legend(loc='upper center', bbox_to_anchor=(0.5, 0),ncol=2)

		clean.errorbar(self.lc.lcs[0].t.loc[good_ix,'MJD'], self.lc.lcs[0].t.loc[good_ix,'uJy'], yerr=self.lc.lcs[0].t.loc[good_ix,'duJy'], fmt='none',ecolor=color,elinewidth=1,c=color)
		clean.scatter(self.lc.lcs[0].t.loc[good_ix,'MJD'], self.lc.lcs[0].t.loc[good_ix,'uJy'], s=50,color=color,marker='o',label='Kept measurements')
		clean.set_title('Kept measurements only')
		clean.axhline(linewidth=1,color='k')
		clean.set_xlabel('MJD')
		clean.set_ylabel('Flux (uJy)')

		self.pdf.savefig(fig, bbox_inches='tight')

	def plot_limcuts(self, limcuts, contam_cut, loss_cut, contam_lim, loss_lim, min_cut, max_cut):
		fig = plt.figure(figsize=(10,6), tight_layout=True)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().spines['top'].set_visible(False)
		plt.title(f'SN {self.lc.tnsname} {self.filt}-band chi-square cut')

		plt.axhline(linewidth=1,color='k')
		plt.xlabel('Chi-square cut')
		plt.ylabel(f'% of baseline measurements')

		plt.axhline(loss_lim,linewidth=1,color='r',linestyle='--',label='Loss limit')
		plt.plot(limcuts.t['PSF Chi-Square Cut'], limcuts.t['Ploss'],ms=5,color='r',marker='o',label='Loss')
		plt.axvline(x=loss_cut,color='r',label='Loss cut')
		plt.axvspan(loss_cut, max_cut, alpha=0.2, color='r')

		plt.axhline(contam_lim,linewidth=1,color='g',linestyle='--',label='Contamination limit')
		plt.plot(limcuts.t['PSF Chi-Square Cut'], limcuts.t['Pcontamination'],ms=5,color='g',marker='o',label='Contamination')
		plt.axvline(x=contam_cut,color='g',label='Contamination cut')
		plt.axvspan(min_cut, contam_cut, alpha=0.2, color='g')
		
		plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')

		self.pdf.savefig(fig, bbox_inches='tight')

	def plot_sim_bumps(self, simparams=None, add2title=None):
		if not self.lc.is_averaged:
			raise RuntimeWarning('ERROR: Light curve to be plotted is not averaged!')

		fig = plt.figure(figsize=(10,6), tight_layout=True)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().spines['top'].set_visible(False)
		plt.axhline(linewidth=1,color='k')
		plt.ylabel('Flux (µJy)')
		plt.xlabel('MJD')
		title = f'SN {self.lc.tnsname} {self.filt}-band averaged flux'
		if not(simparams is None):
			title += f'\nand {len(simparams["sim_peakMJD"].split(","))} simulated pre-SN bump(s) with appmag {simparams["sim_appmag"]:0.2f}'
		if not(add2title is None):
			title += add2title
		plt.title(title)
		plt.axvline(x=self.tchange1, color='magenta', label='ATLAS template change')
		plt.axvline(x=self.tchange2, color='magenta')

		color = 'orange' if self.filt == 'o' else 'cyan'
		good_ix = self.lc.lcs[0].ix_unmasked('Mask',maskval=self.flags['avg_badday'])

		if not(simparams is None):
			plt.errorbar(self.lc.lcs[0].t.loc[good_ix,'MJDbin'], self.lc.lcs[0].t.loc[good_ix,'uJysim'], yerr=self.lc.lcs[0].t.loc[good_ix,'duJy'], zorder=0, fmt='none', ecolor='red', elinewidth=1, c=color)
			plt.scatter(self.lc.lcs[0].t.loc[good_ix,'MJDbin'], self.lc.lcs[0].t.loc[good_ix,'uJysim'], zorder=5, s=45, color='red', marker='o', label=f'simulated bump(s), width={simparams["sim_sigma_plus"]:0.2f} days')
			plt.plot(self.lc.lcs[0].t['MJDbin'], self.lc.lcs[0].t['simLC'], zorder=20, color='red', label='gaussian model(s)')
		plt.errorbar(self.lc.lcs[0].t.loc[good_ix,'MJDbin'], self.lc.lcs[0].t.loc[good_ix,'uJy'], yerr=self.lc.lcs[0].t.loc[good_ix,'duJy'], zorder=10, fmt='none', ecolor=color, elinewidth=1, c=color)
		plt.scatter(self.lc.lcs[0].t.loc[good_ix,'MJDbin'], self.lc.lcs[0].t.loc[good_ix,'uJy'], zorder=15, s=45, color=color, marker='o', label='kept measurements')

		plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)

		xlim_lower = self.args.xlim_lower if not(self.args.xlim_lower is None) else self.xlim_lower
		xlim_upper = self.args.xlim_upper if not(self.args.xlim_upper is None) else self.xlim_upper
		plt.xlim(xlim_lower, xlim_upper)
		ylim_lower = self.args.ylim_lower if not(self.args.ylim_lower is None) else self.ylim_lower
		ylim_upper = self.args.ylim_upper if not(self.args.ylim_upper is None) else self.ylim_upper
		plt.ylim(ylim_lower, ylim_upper)

		self.pdf.savefig(fig)

	def plot_snr(self, simparams=None, add2title=None):
		if not self.lc.is_averaged:
			raise RuntimeWarning('ERROR: Light curve to be plotted is not averaged!')

		fig = plt.figure(figsize=(10,6), tight_layout=True)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().spines['top'].set_visible(False)
		plt.axhline(linewidth=1,color='k')
		plt.ylabel('Signal-to-noise')
		plt.xlabel('MJD')
		title = f'SN {self.lc.tnsname} {self.filt}-band signal-to-noise \nand gaussian weighted rolling sum of signal-to-noise'
		if not(add2title is None):
			title += add2title
		plt.title(title)

		color = 'orange' if self.filt == 'o' else 'cyan'
		good_ix = self.lc.lcs[0].ix_unmasked('Mask',maskval=self.flags['avg_badday'])

		if not(simparams is None):
			plt.scatter(self.lc.lcs[0].t.loc[good_ix,'MJDbin'], self.lc.lcs[0].t.loc[good_ix,'SNRsim'], s=45, color='red', marker='o', zorder=0, label='simulated S/N')
			plt.plot(self.lc.lcs[0].t['MJDbin'], self.lc.lcs[0].t['SNRsimsum'], color='red', zorder=15, label='gaussian weighted rolling sum of simulated S/N')
		plt.scatter(self.lc.lcs[0].t.loc[good_ix,'MJDbin'], self.lc.lcs[0].t.loc[good_ix,'SNR'], s=45, color=color, marker='o', zorder=5, label='S/N')
		plt.plot(self.lc.lcs[0].t['MJDbin'], self.lc.lcs[0].t['SNRsumnorm'], color=color, zorder=10, label='gaussian weighted rolling sum of S/N')
		
		plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)

		self.pdf.savefig(fig)

	def plot_all_snr(self, simparams=None, add2title=None):
		if not self.lc.is_averaged:
			raise RuntimeWarning('ERROR: Light curve to be plotted is not averaged!')

		fig = plt.figure(figsize=(10,6), tight_layout=True)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().spines['top'].set_visible(False)
		plt.axhline(linewidth=1,color='k')
		plt.ylabel('Signal-to-noise')
		plt.xlabel('MJD')
		title = f'SN {self.lc.tnsname} and control light curves {self.filt}-band \ngaussian weighted rolling sum of signal-to-noise'
		if not(add2title is None):
			title += add2title
		plt.title(title)

		for control_index in self.lc.lcs:
			label = None
			if control_index == 0:
				if not(simparams is None):
					color = 'red'
				elif self.filt == 'o':
					color = 'orange'
				else:
					color = 'cyan'

				label = f'SN {self.lc.tnsname}' if simparams is None else f'SN {self.lc.tnsname} with simulated bump(s)'
				zorder = 10
				snrsum_colname = 'SNRsumnorm' if simparams is None else 'SNRsimsum'
			else:
				color = 'blue'
				if control_index == 1:
					label = f'{len(self.lc.lcs)-1} control light curve(s)'
				zorder = 0
				snrsum_colname = 'SNRsumnorm'

			plt.plot(self.lc.lcs[control_index].t['MJDbin'], self.lc.lcs[control_index].t[snrsum_colname], color=color, zorder=zorder, label=label)

		plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)

		self.pdf.savefig(fig)
