#!/usr/bin/env python

import argparse
from typing import List
from download import load_config
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
    color=["green", "blue", "purple", "magenta"]
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

# color scheme
SN_FLUX_COLORS = {"o": "orange", "c": "cyan"}
SN_FLAGGED_FLUX_COLOR = "red"
CONTROL_FLUX_COLOR = "steelblue"

# ATLAS template change dates
TEMPLATE_CHANGE_1_MJD = 58417
TEMPLATE_CHANGE_2_MJD = 58882


class PlotLimits:
    def __init__(self, xlower=None, xupper=None, ylower=None, yupper=None):
        self.xlower = xlower
        self.xupper = xupper
        self.ylower = ylower
        self.yupper = yupper

    def set_lims(self, xlower=None, xupper=None, ylower=None, yupper=None):
        self.xlower = xlower
        self.xupper = xupper
        self.ylower = ylower
        self.yupper = yupper

    def calc_lims(self, lc: LightCurve, indices: List[int] = None) -> List[int]:
        pass

    def get_xlims(self):
        return [self.xlower, self.xupper]

    def get_ylims(self):
        return [self.ylower, self.yupper]

    def __str__(self):
        return f"Plot limits: x-axis [{self.xlower}, {self.xupper}], y-axis [{self.ylower}, {self.yupper}]"


class Plot:
    def __init__(self):
        pass

    def get_lims(self, lc: LightCurve, indices: List[int] = None) -> PlotLimits:
        lims = PlotLimits()
        lims.calc_lims(lc, indices=indices)
        return lims

    def plot_SN(
        self,
        sn: Supernova,
        lims: PlotLimits,
        plot_controls: bool = True,
        plot_template_changes: bool = True,
    ):
        fig, ax1 = plt.subplots(1, constrained_layout=True)
        fig.set_figwidth(7)
        fig.set_figheight(4)

        title = f"SN {sn.tnsname}"
        if sn.num_controls > 0:
            title += f" & control light curves"
        title += f" {sn.filt}-band flux"
        ax1.set_title(title)

        ax1.minorticks_on()
        ax1.tick_params(direction="in", which="both")
        ax1.set_ylabel(r"Flux ($\mu$Jy)")
        ax1.set_xlabel("MJD")
        ax1.axhline(linewidth=1, color="k")

        if plot_controls and sn.num_controls > 0:
            # plot control light curves
            label = f"{sn.num_controls} control light curves"
            for control_index in sn.get_control_indices():
                lc = sn.lcs[control_index]

                plt.errorbar(
                    lc.t["MJD"],
                    lc.t["uJy"],
                    yerr=lc.t[lc.dflux_colname],
                    fmt="none",
                    ecolor=CONTROL_FLUX_COLOR,
                    elinewidth=1.5,
                    capsize=1.2,
                    c=CONTROL_FLUX_COLOR,
                    alpha=0.5,
                    zorder=0,
                )
                plt.scatter(
                    lc.t["MJD"],
                    lc.t["uJy"],
                    s=marker_size,
                    color=CONTROL_FLUX_COLOR,
                    marker="o",
                    alpha=0.5,
                    zorder=0,
                    label=label,
                )

                if not label is None:
                    label = None

        lc = sn.lcs[0]
        preMJD0_ix = lc.get_preMJD0_indices(sn.mjd0)
        postMJD0_ix = lc.get_postMJD0_indices(sn.mjd0)

        # plot pre-MJD0 SN light curve
        plt.errorbar(
            lc.t.loc[preMJD0_ix, "MJD"],
            lc.t.loc[preMJD0_ix, "uJy"],
            yerr=lc.t.loc[preMJD0_ix, lc.dflux_colname],
            fmt="none",
            ecolor="magenta",
            elinewidth=1,
            capsize=1.2,
            c="magenta",
            alpha=0.5,
            zorder=10,
        )
        plt.scatter(
            lc.t.loc[preMJD0_ix, "MJD"],
            lc.t.loc[preMJD0_ix, "uJy"],
            s=marker_size,
            lw=marker_edgewidth,
            color="magenta",
            marker="o",
            alpha=0.5,
            zorder=10,
            label="Pre-MJD0 light curve",
        )

        # plot post-MJD0 SN light curve
        plt.errorbar(
            lc.t.loc[postMJD0_ix, "MJD"],
            lc.t.loc[postMJD0_ix, "uJy"],
            yerr=lc.t.loc[postMJD0_ix, lc.dflux_colname],
            fmt="none",
            ecolor="lime",
            elinewidth=1,
            capsize=1.2,
            c="lime",
            alpha=0.5,
            zorder=10,
        )
        plt.scatter(
            lc.t.loc[postMJD0_ix, "MJD"],
            lc.t.loc[postMJD0_ix, "uJy"],
            s=marker_size,
            lw=marker_edgewidth,
            color="lime",
            marker="o",
            alpha=0.5,
            zorder=10,
            label="Post-MJD0 light curve",
        )

        if plot_template_changes:
            pass
            # TODO

        ax1.set_xlim(lims.xlower, lims.xupper)
        ax1.set_ylim(lims.ylower, lims.yupper)
        ax1.legend(loc="upper right", facecolor="white", framealpha=1.0).set_zorder(100)

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


class PlotPdf(Plot):
    def __init__(self, output_dir, tnsname, filt="o"):
        Plot.__init__(self)
        filename = f"{output_dir}/{tnsname}_plots.{filt}.pdf"
        self.pdf = PdfPages(filename)

    def save_pdf(self):
        print("\nSaving PDF of plots...\n")
        self.pdf.close()

    def save_fig(self, figure: Figure):
        self.pdf.savefig(figure)


# define command line arguments
def define_args(parser=None, usage=None, conflict_handler="resolve"):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

    parser.add_argument(
        "tnsnames", nargs="+", help="TNS names of the transients to clean"
    )
    parser.add_argument(
        "--config_file",
        default="config.ini",
        type=str,
        help="file name of .ini file with settings for this class",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        default=False,
        action="store_true",
        help="overwrite existing file with same file name",
    )
    parser.add_argument(
        "--filters",
        type=str,
        default=None,
        help="comma-separated list of filters to plot",
    )
    parser.add_argument(
        "--num_controls",
        type=int,
        default=None,
        help="number of control light curves to load and plot",
    )

    # possible cuts
    parser.add_argument(
        "-t",
        "--template_correction",
        default=False,
        action="store_true",
        help="plot ATLAS template change correction",
    )
    parser.add_argument(
        "-e",
        "--uncert_est",
        default=False,
        action="store_true",
        help="plot true uncertainty estimation",
    )
    parser.add_argument(
        "-u",
        "--uncert_cut",
        default=False,
        action="store_true",
        help="plot uncertainty cut",
    )
    parser.add_argument(
        "-x",
        "--x2_cut",
        default=False,
        action="store_true",
        help="plot chi-square cut",
    )
    parser.add_argument(
        "-c",
        "--controls_cut",
        default=False,
        action="store_true",
        help="plot control light curve cut",
    )
    parser.add_argument(
        "-g",
        "--averaging",
        default=False,
        action="store_true",
        help="plot averaged light curves and bad day cut",
    )
    parser.add_argument(
        "--custom_cuts",
        default=False,
        action="store_true",
        help="scan config file for custom cuts and plot them",
    )

    return parser


if __name__ == "__main__":
    args = define_args().parse_args()
    config = load_config(args.config_file)
