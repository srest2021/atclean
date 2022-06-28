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

	def set(self, lc, filt):
		self.lc = lc 
		self.filt = filt

	def save(self):
		self.pdf.close()

	def plot_lc(self):

		self.pdf.savefig(fig)

	def plot_cut_lc(self, fig, flags, is_averaged=False):

		return fig

	def plot_averaged_lc(self):
		fig = plt.figure()
		fig = self.plot_cut_lc(fig, flags= , is_averaged=True)
		self.pdf.savefig(fig)

	def plot_chisquare_cut(self):

		self.pdf.savefig(fig)

	def plot_uncertainty_cut(self):

		self.pdf.savefig(fig)

	def plot_controls_cut(self):

		self.pdf.savefig(fig)
