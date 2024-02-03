from pdastro import AorB
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

import warnings
warnings.simplefilter('error', RuntimeWarning)
warnings.filterwarnings("ignore")

# plotting styles
plt.rc('axes', titlesize = 17)
plt.rc('xtick', labelsize = 12)
plt.rc('ytick', labelsize = 12)
plt.rc('legend', fontsize = 10)
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=['red', 'orange', 'green', 'blue', 'purple', 'magenta'])
matplotlib.rcParams['xtick.major.size'] = 6
matplotlib.rcParams['xtick.major.width'] = 1
matplotlib.rcParams['xtick.minor.size'] = 3
matplotlib.rcParams['xtick.minor.width'] = 1
matplotlib.rcParams['ytick.major.size'] = 6
matplotlib.rcParams['ytick.major.width'] = 1
matplotlib.rcParams['ytick.minor.size'] = 3
matplotlib.rcParams['ytick.minor.width'] = 1
matplotlib.rcParams['axes.linewidth'] = 1
marker_size = 30
marker_edgewidth = 1.5
sn_flagged_flux = 'red' 
ctrl_flux = 'steelblue'

class PlotAtlasLightCurve():
	def __init__(self, lc, flags, add2filename=None):
		self.lc = lc
		self.flags = flags

		output_dir = self.lc.cfg['output_dir']
		if add2filename is None:
			self.pdf = PdfPages(f'{output_dir}/{lc.tnsname}/{lc.tnsname}_{lc.filt}_plots.pdf')
		else:
			self.pdf = PdfPages(f'{output_dir}/{lc.tnsname}/{lc.tnsname}_{lc.filt}_plots_{add2filename}.pdf')

		if not 'plot_params' in self.lc.cfg:
			raise RuntimeError('ERROR: attempting to plot, but no plot params!')
		manual_limits = list(self.lc.cfg['plot_params'].values()) 
		self.limits = self.set_limits(self.lc, manual_limits)

		# ATLAS template change dates
		self.t1 = 58417
		self.t2 = 58882

	# save PDF of current plots
	def save(self):
		print('\nSaving PDF of plots...\n')
		self.pdf.close()

	def set_limits(self, lc, manual_limits, indices=None):
		limits = manual_limits
		if limits[0] is None:
			limits[0] = lc.lcs[0].t['MJD'].min() * 0.999
		if limits[1] is None:
			limits[1] = lc.lcs[0].t['MJD'].max() * 1.001
		
		if indices is None:
			indices = lc.get_ix()
		# exclude measurements with duJy > 160
		good_ix = lc.lcs[0].ix_inrange(colnames='duJy', uplim=160, indices=indices)
		# get 5% of abs(max flux - min flux)
		flux_min = lc.lcs[0].t.loc[good_ix, 'uJy'].min()
		flux_max = lc.lcs[0].t.loc[good_ix, 'uJy'].max()
		diff = abs(flux_max - flux_min)
		offset = 0.05 * diff

		if limits[2] is None: 
			limits[2] = flux_min - offset
		if limits[3] is None:
			limits[3] = flux_max + offset

		return limits
		
	def plot_lcs(self, add2title=None, plot_templates=True, plot_controls=True):
		sn_flux = 'orange' if self.lc.filt == 'o' else 'cyan'
		fig, ax1 = plt.subplots(1, constrained_layout=True)
		fig.set_figwidth(7)
		fig.set_figheight(4)

		title = f'SN {self.lc.tnsname} & control light curves {self.lc.filt}-band flux'
		if not(add2title is None):
			title += add2title
		ax1.set_title(title)

		ax1.minorticks_on()
		ax1.tick_params(direction='in', which='both')
		ax1.set_ylabel(r'Flux ($\mu$Jy)')
		ax1.set_xlabel('MJD')
		ax1.axhline(linewidth=1, color='k')

		ax1.set_xlim(self.limits[0],self.limits[1])
		ax1.set_ylim(self.limits[2],self.limits[3])

		if plot_templates:
			ax1.axvline(x=self.t1, color='k', linestyle='dotted', label='ATLAS template change', zorder=100)
			ax1.axvline(x=self.t2, color='k', linestyle='dotted', zorder=100)

		preSN_ix = self.lc.get_pre_SN_ix()
		postSN_ix = self.lc.get_post_SN_ix()

		if plot_controls:
			for control_index in range(1, self.lc.num_controls+1):
				ax1.errorbar(self.lc.lcs[control_index].t['MJD'], self.lc.lcs[control_index].t['uJy'], yerr=self.lc.lcs[control_index].t[self.lc.dflux_colnames[control_index]], fmt='none', ecolor=ctrl_flux, elinewidth=1.5, capsize=1.2, c=ctrl_flux, alpha=0.5, zorder=0)
				if control_index == 1:
					ax1.scatter(self.lc.lcs[control_index].t['MJD'], self.lc.lcs[control_index].t['uJy'], s=marker_size, color=ctrl_flux, marker='o', alpha=0.5, zorder=0, label=f'{self.lc.num_controls} control light curves')
				else:
					ax1.scatter(self.lc.lcs[control_index].t['MJD'], self.lc.lcs[control_index].t['uJy'], s=marker_size, color=ctrl_flux, marker='o', alpha=0.5, zorder=0)


		ax1.errorbar(self.lc.lcs[0].t.loc[preSN_ix,'MJD'], self.lc.lcs[0].t.loc[preSN_ix,'uJy'], yerr=self.lc.lcs[0].t.loc[preSN_ix,self.lc.dflux_colnames[0]], fmt='none', ecolor=sn_flux, elinewidth=1, capsize=1.2, c=sn_flux, alpha=0.5, zorder=10)
		ax1.scatter(self.lc.lcs[0].t.loc[preSN_ix,'MJD'], self.lc.lcs[0].t.loc[preSN_ix,'uJy'], s=marker_size, lw=marker_edgewidth, color=sn_flux, marker='o', alpha=0.5, zorder=10, label='Pre-SN')

		ax1.errorbar(self.lc.lcs[0].t.loc[postSN_ix,'MJD'], self.lc.lcs[0].t.loc[postSN_ix,'uJy'], yerr=self.lc.lcs[0].t.loc[postSN_ix,self.lc.dflux_colnames[0]], fmt='none', ecolor='red', elinewidth=1, capsize=1.2, c=sn_flux, alpha=0.5, zorder=10)
		ax1.scatter(self.lc.lcs[0].t.loc[postSN_ix,'MJD'], self.lc.lcs[0].t.loc[postSN_ix,'uJy'], s=marker_size, lw=marker_edgewidth, color='red', marker='o', alpha=0.5, zorder=10, label='During and after SN')

		ax1.legend(loc='upper right', facecolor='white', framealpha=1.0).set_zorder(100)

		self.pdf.savefig(fig)

	def plot_cut_lc(self, lc, title, flag):
		sn_flux = 'orange' if lc.filt == 'o' else 'cyan'
		fig, (ax2, ax1) = plt.subplots(2, constrained_layout=True)
		fig.set_figwidth(7)
		fig.set_figheight(5)

		fig.suptitle(f'{title} (flag {hex(flag)})')

		ax1.minorticks_on()
		ax1.tick_params(direction='in', which='both')
		ax2.get_xaxis().set_ticks([])
		ax1.set_ylabel(r'Flux ($\mu$Jy)')
		ax1.axhline(linewidth=1, color='k')

		ax2.minorticks_on()
		ax2.tick_params(direction='in', which='both')
		ax2.set_ylabel(r'Flux ($\mu$Jy)')
		ax1.set_xlabel('MJD')
		ax2.axhline(linewidth=1, color='k')

		ax1.set_xlim(self.limits[0],self.limits[1])
		ax1.set_ylim(self.limits[2],self.limits[3])
		ax2.set_xlim(self.limits[0],self.limits[1])
		ax2.set_ylim(self.limits[2],self.limits[3])

		good_ix = lc.lcs[0].ix_unmasked('Mask', maskval=flag)
		bad_ix = lc.lcs[0].ix_masked('Mask', maskval=flag)

		ax1.errorbar(lc.lcs[0].t.loc[good_ix,'MJD'], lc.lcs[0].t.loc[good_ix,'uJy'], yerr=lc.lcs[0].t.loc[good_ix,lc.dflux_colnames[0]], fmt='none', ecolor=sn_flux, elinewidth=1, capsize=1.2, c=sn_flux, alpha=0.5)
		ax1.scatter(lc.lcs[0].t.loc[good_ix,'MJD'], lc.lcs[0].t.loc[good_ix,'uJy'], s=marker_size, lw=marker_edgewidth, color=sn_flux, marker='o', alpha=0.5, label='Kept measurements')

		ax2.errorbar(lc.lcs[0].t.loc[good_ix,'MJD'], lc.lcs[0].t.loc[good_ix,'uJy'], yerr=lc.lcs[0].t.loc[good_ix,lc.dflux_colnames[0]], fmt='none', ecolor=sn_flux, elinewidth=1, capsize=1.2, c=sn_flux, alpha=0.5, zorder=5)
		ax2.scatter(lc.lcs[0].t.loc[good_ix,'MJD'], lc.lcs[0].t.loc[good_ix,'uJy'], s=marker_size, lw=marker_edgewidth, color=sn_flux, marker='o', alpha=0.5, label='Kept measurements', zorder=5)

		ax2.errorbar(lc.lcs[0].t.loc[bad_ix,'MJD'], lc.lcs[0].t.loc[bad_ix,'uJy'], yerr=lc.lcs[0].t.loc[bad_ix,lc.dflux_colnames[0]], fmt='none', ecolor=sn_flagged_flux, elinewidth=1, capsize=1.2, c=sn_flagged_flux, alpha=0.5, zorder=10)
		ax2.scatter(lc.lcs[0].t.loc[bad_ix,'MJD'], lc.lcs[0].t.loc[bad_ix,'uJy'], s=marker_size, lw=marker_edgewidth, color=sn_flagged_flux, facecolors='none', edgecolors=sn_flagged_flux, marker='o', alpha=0.5, label='Cut measurements', zorder=10)

		ax1.legend(loc='upper right', facecolor='white', framealpha=1.0).set_zorder(100)
		ax2.legend(loc='upper right', facecolor='white', framealpha=1.0).set_zorder(100)

		self.pdf.savefig(fig)

	def plot_limcuts(self, limcuts, cut, cut_start, cut_stop, use_preSN_lc=False):
		loss_color = 'darkmagenta'
		contam_color = 'teal'

		fig, ax1 = plt.subplots(1, constrained_layout=True)
		fig.set_figwidth(5.5)
		fig.set_figheight(3)

		ax1.set_title(f'SN {self.lc.tnsname} {self.lc.filt}-band chi-square cut')

		ax1.minorticks_on()
		ax1.tick_params(direction='in', which='both')
		if use_preSN_lc:
			ax1.set_ylabel(f'% pre-SN light curve measurements')
		else:
			ax1.set_ylabel(f'% control light curve measurements')
		ax1.set_xlabel('Chi-square cut')
		ax1.axhline(linewidth=1, color='k')

		ax1.plot(limcuts['PSF Chi-Square Cut'].values, limcuts['Ploss'].values, ms=3.5, color=loss_color, marker='o', label='Loss')
		ax1.plot(limcuts['PSF Chi-Square Cut'].values, limcuts['Pcontamination'].values, ms=3.5, color=contam_color, marker='o', label='Contamination')
		ax1.axvline(x=cut, color='k', linestyle='dashed', label='Selected cut')

		ax1.set_xlim(cut_start-1,cut_stop+1)
		ax1.set_ylim(0, max(max(limcuts['Ploss'].values), max(limcuts['Pcontamination'].values))*1.1)

		ax1.legend(facecolor='white', framealpha=1, bbox_to_anchor=(1.02, 1), loc='upper left')

		self.pdf.savefig(fig)

	def plot_uncert_cut(self, lc):
		self.plot_cut_lc(lc, 'Uncertainty cut', self.flags['uncertainty'])
		
	def plot_x2_cut(self, lc, limcuts, cut, cut_start, cut_stop, use_preSN_lc=False):
		self.plot_limcuts(limcuts, cut, cut_start, cut_stop, use_preSN_lc=use_preSN_lc)
		self.plot_cut_lc(lc, 'Chi-square cut', self.flags['chisquare'])

	def plot_controls_cut(self, lc):
		self.plot_cut_lc(lc, 'Control light curve cut', self.flags['controls_bad'])

	def plot_all_cuts(self, lc):
		self.plot_cut_lc(lc, 'Uncertainty, chi-square, and control light curve cuts', 
						 self.flags['uncertainty']|self.flags['chisquare']|self.flags['controls_bad'])

	def plot_badday_cut(self, avglc):
		self.plot_cut_lc(avglc, 'Averaged light curve', self.flags['avg_badday'])

	def plot_uncert_est(self, lc):
		if not 'duJy_new' in lc.lcs[0].t.columns:
			raise RuntimeError('ERROR: attempting to plot uncertainty estimations, but no "duJy_new" column!')

		sn_flux = 'orange' if lc.filt == 'o' else 'cyan'
		fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True)
		fig.set_figwidth(7)
		fig.set_figheight(5)

		ax1.set_title(f'SN {lc.tnsname} {lc.filt}-band flux\nbefore true uncertainties estimation')
		ax1.minorticks_on()
		ax1.tick_params(direction='in', which='both')
		ax1.get_xaxis().set_ticks([])
		ax1.set_ylabel(r'Flux ($\mu$Jy)')
		ax1.axhline(linewidth=1, color='k')

		ax2.set_title(f'after true uncertainties estimation')
		ax2.minorticks_on()
		ax2.tick_params(direction='in', which='both')
		ax2.set_ylabel(r'Flux ($\mu$Jy)')
		ax2.set_xlabel('MJD')
		ax2.axhline(linewidth=1, color='k')

		ax1.set_xlim(self.limits[0],self.limits[1])
		ax1.set_ylim(self.limits[2],self.limits[3])

		ax1.errorbar(lc.lcs[0].t['MJD'], lc.lcs[0].t['uJy'], yerr=lc.lcs[0].t['duJy'], fmt='none', ecolor=sn_flux, elinewidth=1, capsize=1.2, c=sn_flux, alpha=0.5)
		ax1.scatter(lc.lcs[0].t['MJD'], lc.lcs[0].t['uJy'], s=marker_size, lw=marker_edgewidth, color=sn_flux, marker='o', alpha=0.5)

		ax2.errorbar(lc.lcs[0].t['MJD'], lc.lcs[0].t['uJy'], yerr=lc.lcs[0].t['duJy_new'], fmt='none', ecolor=sn_flux, elinewidth=1, capsize=1.2, c=sn_flux, alpha=0.5)
		ax2.scatter(lc.lcs[0].t['MJD'], lc.lcs[0].t['uJy'], s=marker_size, lw=marker_edgewidth, color=sn_flux, marker='o', alpha=0.5)

		self.pdf.savefig(fig)
	
	def plot_template_correction(self, lc, title=None):
		sn_flux = 'orange' if lc.filt == 'o' else 'cyan'
		colors = ['salmon', 'sandybrown', 'darkseagreen']
		t1, t2 = 58417, 58882

		region1_ix = lc.lcs[0].ix_inrange('MJD', uplim=t1)
		region2_ix = lc.lcs[0].ix_inrange('MJD', lowlim=t1, uplim=t2)
		region3_ix = lc.lcs[0].ix_inrange('MJD', lowlim=t2)

		region1_mean = lc._get_mean(region1_ix[-40:]) # last 40 measurements before t1
		region2a_mean = lc._get_mean(region2_ix[:40]) # first 40 measurements after t1
		region2b_mean = lc._get_mean(region2_ix[-40:]) # last 40 measurements before t2
		region3_mean = lc._get_mean(region3_ix[:40]) # first 40 measurements after t2

		gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], hspace=0.35, wspace=0.4)
		fig = plt.figure(constrained_layout=True)
		fig.set_figwidth(6)
		fig.set_figheight(6)
		fig.tight_layout()

		ax1 = plt.subplot(gs[0, :])
		if not title is None:
			ax1.set_title(title)
		ax1.axvline(x=t1, color='k', linestyle='dotted', label='ATLAS template change', zorder=30)
		ax1.axvline(x=t2, color='k', linestyle='dotted', zorder=30)
		ax1.axhline(color='k',zorder=0)
		limits = self.set_limits(lc, [None]*4)
		ax1.set_xlim(limits[0], limits[1])
		ax1.set_ylim(limits[2], limits[3])

		ax1.errorbar(lc.lcs[0].t.loc[region1_ix,'MJD'], lc.lcs[0].t.loc[region1_ix,'uJy'], yerr=lc.lcs[0].t.loc[region1_ix,lc.dflux_colnames[0]], fmt='none', ecolor=colors[0], elinewidth=1, capsize=1.2, c=sn_flux, alpha=0.5, zorder=10)
		ax1.scatter(lc.lcs[0].t.loc[region1_ix,'MJD'], lc.lcs[0].t.loc[region1_ix,'uJy'], s=marker_size, lw=marker_edgewidth, color=colors[0], marker='o', alpha=0.5, zorder=10, label='Region 1 flux')
		ax1.errorbar(lc.lcs[0].t.loc[region2_ix,'MJD'], lc.lcs[0].t.loc[region2_ix,'uJy'], yerr=lc.lcs[0].t.loc[region2_ix,lc.dflux_colnames[0]], fmt='none', ecolor=colors[1], elinewidth=1, capsize=1.2, c=sn_flux, alpha=0.5, zorder=10)
		ax1.scatter(lc.lcs[0].t.loc[region2_ix,'MJD'], lc.lcs[0].t.loc[region2_ix,'uJy'], s=marker_size, lw=marker_edgewidth, color=colors[1], marker='o', alpha=0.5, zorder=10, label='Region 2 flux')
		ax1.errorbar(lc.lcs[0].t.loc[region3_ix,'MJD'], lc.lcs[0].t.loc[region3_ix,'uJy'], yerr=lc.lcs[0].t.loc[region3_ix,lc.dflux_colnames[0]], fmt='none', ecolor=colors[2], elinewidth=1, capsize=1.2, c=sn_flux, alpha=0.5, zorder=10)
		ax1.scatter(lc.lcs[0].t.loc[region3_ix,'MJD'], lc.lcs[0].t.loc[region3_ix,'uJy'], s=marker_size, lw=marker_edgewidth, color=colors[2], marker='o', alpha=0.5, zorder=10, label='Region 2 flux')

		ax1.legend(facecolor='white', framealpha=1, loc='upper left').set_zorder(100)#,  bbox_to_anchor=(1, 1))

		ax2 = plt.subplot(gs[1, 0])
		ax2.set_title('First template change', fontsize=12)
		ax2.axvline(x=t1, color='k', linestyle='dotted', zorder=100)
		ax2.axhline(color='k',zorder=0)
		ax2.set_xlim(lc.lcs[0].t.loc[region1_ix[-40:][0], 'MJD'], lc.lcs[0].t.loc[region2_ix[:40][-1], 'MJD'])

		ax2.errorbar(lc.lcs[0].t.loc[region1_ix,'MJD'], lc.lcs[0].t.loc[region1_ix,'uJy'], yerr=lc.lcs[0].t.loc[region1_ix,lc.dflux_colnames[0]], fmt='none', ecolor=colors[0], elinewidth=1, capsize=1.2, c=sn_flux, alpha=0.5, zorder=10)
		ax2.scatter(lc.lcs[0].t.loc[region1_ix,'MJD'], lc.lcs[0].t.loc[region1_ix,'uJy'], s=marker_size, lw=marker_edgewidth, color=colors[0], marker='o', alpha=0.5, zorder=10)
		ax2.errorbar(lc.lcs[0].t.loc[region2_ix,'MJD'], lc.lcs[0].t.loc[region2_ix,'uJy'], yerr=lc.lcs[0].t.loc[region2_ix,lc.dflux_colnames[0]], fmt='none', ecolor=colors[1], elinewidth=1, capsize=1.2, c=sn_flux, alpha=0.5, zorder=10)
		ax2.scatter(lc.lcs[0].t.loc[region2_ix,'MJD'], lc.lcs[0].t.loc[region2_ix,'uJy'], s=marker_size, lw=marker_edgewidth, color=colors[1], marker='o', alpha=0.5, zorder=10)

		ax2.axhline(y=region1_mean, color=colors[0], linestyle='dashed', label='Region 1 mean',zorder=20)
		ax2.axhline(y=region2a_mean, color=colors[1], linestyle='dashed', label='Region 2 mean',zorder=20)

		limits = self.set_limits(lc, [None]*4, indices=AorB(region1_ix[-40:], region2_ix[:40]))
		ax2.set_ylim(limits[2], limits[3])
		ax2.legend(facecolor='white', framealpha=1).set_zorder(100)

		ax3 = plt.subplot(gs[1, 1]) 
		ax3.set_title('Second template change', fontsize=12)
		ax3.axvline(x=t2, color='k', linestyle='dotted', zorder=30)
		ax3.axhline(color='k',zorder=0)
		ax3.set_xlim(lc.lcs[0].t.loc[region2_ix[-40:][0], 'MJD'], lc.lcs[0].t.loc[region3_ix[:40][-1], 'MJD'])
		ax3.set_ylim(limits[2], limits[3])

		ax3.errorbar(lc.lcs[0].t.loc[region2_ix,'MJD'], lc.lcs[0].t.loc[region2_ix,'uJy'], yerr=lc.lcs[0].t.loc[region2_ix,lc.dflux_colnames[0]], fmt='none', ecolor=colors[1], elinewidth=1, capsize=1.2, c=sn_flux, alpha=0.5, zorder=10)
		ax3.scatter(lc.lcs[0].t.loc[region2_ix,'MJD'], lc.lcs[0].t.loc[region2_ix,'uJy'], s=marker_size, lw=marker_edgewidth, color=colors[1], marker='o', alpha=0.5, zorder=10)
		ax3.errorbar(lc.lcs[0].t.loc[region3_ix,'MJD'], lc.lcs[0].t.loc[region3_ix,'uJy'], yerr=lc.lcs[0].t.loc[region3_ix,lc.dflux_colnames[0]], fmt='none', ecolor=colors[2], elinewidth=1, capsize=1.2, c=sn_flux, alpha=0.5, zorder=10)
		ax3.scatter(lc.lcs[0].t.loc[region3_ix,'MJD'], lc.lcs[0].t.loc[region3_ix,'uJy'], s=marker_size, lw=marker_edgewidth, color=colors[2], marker='o', alpha=0.5, zorder=10)

		ax3.axhline(y=region2b_mean, color=colors[1], linestyle='dashed', label='Region 2 mean',zorder=20)
		ax3.axhline(y=region3_mean, color=colors[2], linestyle='dashed', label='Region 3 mean',zorder=20)

		limits = self.set_limits(lc, [None]*4, indices=AorB(region2_ix[-40:], region3_ix[:40]))
		ax3.set_ylim(limits[2], limits[3])
		ax3.legend(facecolor='white', framealpha=1).set_zorder(100)

		for ax in (ax1, ax2, ax3):
			ax.minorticks_on()
			ax.tick_params(direction='in', which='both')
			ax.set_xlabel('MJD')
			ax.set_ylabel('uJy')

		self.pdf.savefig(fig)