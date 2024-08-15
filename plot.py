#!/usr/bin/env python

from typing import List
from pdastro import AorB
import matplotlib
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
from lightcurve import LightCurve, AveragedLightCurve, Supernova, AveragedSupernova
from clean import Cut, LimCutsTable

# plotting styles
plt.rc("axes", titlesize=17)
plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=12)
plt.rc("legend", fontsize=10)
plt.rcParams["font.size"] = 12
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["axes.prop_cycle"] = matplotlib.cycler(
    color=["red", "orange", "green", "blue", "purple", "magenta"]
)
matplotlib.rcParams["xtick.major.size"] = 6
matplotlib.rcParams["xtick.major.width"] = 1
matplotlib.rcParams["xtick.minor.size"] = 3
matplotlib.rcParams["xtick.minor.width"] = 1
matplotlib.rcParams["ytick.major.size"] = 6
matplotlib.rcParams["ytick.major.width"] = 1
matplotlib.rcParams["ytick.minor.size"] = 3
matplotlib.rcParams["ytick.minor.width"] = 1
matplotlib.rcParams["axes.linewidth"] = 1
marker_size = 30
marker_edgewidth = 1.5
sn_flagged_flux = "red"
ctrl_flux = "steelblue"

# ATLAS template change dates
TEMPLATE_CHANGE_1_MJD = 58417
TEMPLATE_CHANGE_2_MJD = 58882


class PlotPdf:
    def __init__(self, output_dir, tnsname, filt="o"):
        filename = f"{output_dir}/{tnsname}_plots.{filt}.pdf"
        self.pdf = PdfPages(filename)

    def save_pdf(self):
        print("\nSaving PDF of plots...\n")
        self.pdf.close()

    def save_fig(self, figure: Figure):
        self.pdf.savefig(figure)

    def get_lims(self, lc: LightCurve, indices: List[int] = None) -> List[int]:
        pass

    def plot_SN(self, sn: Supernova):
        pass

    def plot_averaged_SN(self, sn: AveragedSupernova):
        pass

    def plot_cut(self, lc: LightCurve, cut: Cut, title: str | None = None):
        pass

    def plot_limcuts(
        self,
        limcuts: LimCutsTable,
        cut: Cut,
        cut_start: int,
        cut_stop: int,
        use_preSN_lc=False,
    ):
        pass

    def plot_uncert_est(self, lc: LightCurve):
        pass

    def plot_template_correction(self, lc: LightCurve):
        pass
