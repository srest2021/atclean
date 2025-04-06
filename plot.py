#!/usr/bin/env python

import os
from typing import List, Optional
import matplotlib
from matplotlib import gridspec
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from lightcurve import (
    LimCutsTable,
    LightCurve,
    Supernova,
    AveragedSupernova,
)
from utils import TEMPLATE_CHANGE_1_MJD, TEMPLATE_CHANGE_2_MJD, ChiSquareCut, PlotLimits

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
SN_FLUX_COLORS = {
    # ATLAS
    "o": "orange",  # Orange filter
    "c": "cyan",  # Cyan filter
    # Rubin
    "u": "purple",  # Ultraviolet (u-band)
    "g": "green",  # Green (g-band)
    "r": "salmon",  # Red (r-band)
    "i": "indigo",  # Near-infrared (i-band)
    "z": "brown",  # Deep red (z-band)
    "y": "darkred",  # Near-infrared (y-band)
    # TESS
    "tess": "pink",  # TESS uses a single wide bandpass (red-sensitive)
}
SN_FLAGGED_FLUX_COLOR = "red"
CONTROL_FLUX_COLOR = "steelblue"
SELECT_CONTROL_FLUX_COLOR = "forestgreen"
FACE_COLOR = "whitesmoke"
SN_FOM_COLOR = "deeppink"
CONTROL_FOM_COLOR = "cornflowerblue"
SELECT_CONTROL_FOM_COLOR = "mediumblue"
colors = [
    "indianred",
    "salmon",
    "sandybrown",
    "gold",
    "yellowgreen",
    "mediumseagreen",
    "turquoise",
    "lightskyblue",
    "plum",
    "palevioletred",
]
plt.rcParams["axes.prop_cycle"] = matplotlib.cycler(color=colors)

# line styles
SIM_BUMP_LS = "dashed"
FOM_LIMIT_LS = "dotted"


class Plot:
    def __init__(self, output_dir: str = None):
        self.output_dir = output_dir

    def save_plot(self, filename, **kwargs):
        filename = f"{self.output_dir}/{filename}.png"
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        print(f"Saving plot: {filename}")
        plt.savefig(filename, dpi=200, **kwargs)

    def get_lims(
        self,
        sn: Supernova = None,
        control_index: int = 0,
        flag: Optional[int] = None,
        custom_lims: Optional[PlotLimits] = None,
        indices: List[int] = None,
        pre_sn: bool = False,
    ) -> PlotLimits:
        lims = PlotLimits()

        if custom_lims is not None and custom_lims.get_xlims() is not None:
            lims.set_xlims(custom_lims.get_xlims())
        else:
            lims.set_xlims(
                sn.lcs[control_index].get_xlims(
                    mjd0=sn.mjd0 if pre_sn else None,
                )
            )

        if custom_lims is not None and custom_lims.get_ylims() is not None:
            lims.set_ylims(custom_lims.get_ylims())
        else:
            lims.set_ylims(sn.lcs[control_index].get_ylims(indices=indices, flag=flag))

        return lims

    def plot_SN(
        self,
        sn: Supernova,
        custom_lims: Optional[PlotLimits] = None,
        plot_controls: bool = True,
        plot_template_changes: bool = True,
        save: bool = False,
        filename: str = "original",
    ):
        fig, ax1 = plt.subplots(1, constrained_layout=True)
        fig.set_figwidth(7)
        fig.set_figheight(4)

        title = f"SN {sn.tnsname}"
        if plot_controls and sn.num_controls > 0:
            title += f" & control light curves"
        title += f" {sn.filt}-band flux"
        ax1.set_title(title)

        ax1.minorticks_on()
        ax1.tick_params(direction="in", which="both")
        ax1.set_ylabel(r"Flux ($\mu$Jy)")
        ax1.set_xlabel(sn.colnames.mjd)
        ax1.axhline(linewidth=1, color="k")

        if plot_controls and sn.num_controls > 0:
            # plot control light curves
            label = f"{sn.num_controls} control light curves"
            for control_index in sn.get_control_lc_indices():
                lc = sn.lcs[control_index]

                plt.errorbar(
                    lc.t[lc.colnames.mjd],
                    lc.t[lc.colnames.flux],
                    yerr=lc.t[lc.colnames.dflux],
                    fmt="none",
                    ecolor=CONTROL_FLUX_COLOR,
                    elinewidth=1.5,
                    capsize=1.2,
                    c=CONTROL_FLUX_COLOR,
                    alpha=0.5,
                    zorder=0,
                )
                plt.scatter(
                    lc.t[lc.colnames.mjd],
                    lc.t[lc.colnames.flux],
                    s=marker_size,
                    color=CONTROL_FLUX_COLOR,
                    marker="o",
                    alpha=0.5,
                    zorder=0,
                    label=label,
                )

                if not label is None:
                    label = None

        sn_lc = sn.lcs[0]
        preMJD0_ix = sn_lc.get_preMJD0_indices(sn.mjd0)
        postMJD0_ix = sn_lc.get_postMJD0_indices(sn.mjd0)

        if sn_lc.can_plot(preMJD0_ix):
            # plot pre-MJD0 SN light curve
            plt.errorbar(
                sn_lc.t.loc[preMJD0_ix, sn_lc.colnames.mjd],
                sn_lc.t.loc[preMJD0_ix, sn_lc.colnames.flux],
                yerr=sn_lc.t.loc[preMJD0_ix, sn_lc.colnames.dflux],
                fmt="none",
                ecolor="magenta",
                elinewidth=1,
                capsize=1.2,
                c="magenta",
                alpha=0.5,
                zorder=10,
            )
            plt.scatter(
                sn_lc.t.loc[preMJD0_ix, sn_lc.colnames.mjd],
                sn_lc.t.loc[preMJD0_ix, sn_lc.colnames.flux],
                s=marker_size,
                lw=marker_edgewidth,
                color="magenta",
                marker="o",
                alpha=0.5,
                zorder=10,
                label="Pre-MJD0 light curve",
            )

        if sn_lc.can_plot(postMJD0_ix):
            # plot post-MJD0 SN light curve
            plt.errorbar(
                sn_lc.t.loc[postMJD0_ix, sn_lc.colnames.mjd],
                sn_lc.t.loc[postMJD0_ix, sn_lc.colnames.flux],
                yerr=sn_lc.t.loc[postMJD0_ix, sn_lc.colnames.dflux],
                fmt="none",
                ecolor="lime",
                elinewidth=1,
                capsize=1.2,
                c="lime",
                alpha=0.5,
                zorder=10,
            )
            plt.scatter(
                sn_lc.t.loc[postMJD0_ix, sn_lc.colnames.mjd],
                sn_lc.t.loc[postMJD0_ix, sn_lc.colnames.flux],
                s=marker_size,
                lw=marker_edgewidth,
                color="lime",
                marker="o",
                alpha=0.5,
                zorder=10,
                label="Post-MJD0 light curve",
            )

        if plot_template_changes:
            ax1.axvline(
                x=TEMPLATE_CHANGE_1_MJD,
                color="k",
                linestyle="dotted",
                label="ATLAS template change",
                zorder=100,
            )
            ax1.axvline(
                x=TEMPLATE_CHANGE_2_MJD, color="k", linestyle="dotted", zorder=100
            )

        lims = self.get_lims(sn=sn, custom_lims=custom_lims)
        if lims.get_xlims() is not None:
            ax1.set_xlim(lims.get_xlims())
        if lims.get_ylims() is not None:
            ax1.set_ylim(lims.get_ylims())

        ax1.legend(loc="upper right", facecolor="white", framealpha=1.0).set_zorder(100)

        if save:
            self.save_plot(filename)

        return fig

    def plot_cut(
        self,
        sn: Supernova,
        flag: int,
        control_index: int = 0,
        custom_lims: Optional[PlotLimits] = None,
        title: str | None = None,
        save_filename: str = None,
    ):
        fig, (ax2, ax1) = plt.subplots(2, constrained_layout=True)
        fig.set_figwidth(7)
        fig.set_figheight(5)

        if not title:
            title = "Cut"
        fig.suptitle(f"{title} (flag {hex(flag)})")

        ax1.minorticks_on()
        ax1.tick_params(direction="in", which="both")
        ax2.get_xaxis().set_ticks([])
        ax1.set_ylabel(r"Flux ($\mu$Jy)")
        ax1.axhline(linewidth=1, color="k")

        ax2.minorticks_on()
        ax2.tick_params(direction="in", which="both")
        ax2.set_ylabel(r"Flux ($\mu$Jy)")
        ax1.set_xlabel(sn.lcs[control_index].colnames.mjd)
        ax2.axhline(linewidth=1, color="k")

        good_ix = sn.lcs[control_index].get_good_indices(flag)
        bad_ix = sn.lcs[control_index].get_bad_indices(flag)

        if sn.lcs[control_index].can_plot(good_ix):
            ax1.errorbar(
                sn.lcs[control_index].t.loc[
                    good_ix, sn.lcs[control_index].colnames.mjd
                ],
                sn.lcs[control_index].t.loc[
                    good_ix, sn.lcs[control_index].colnames.flux
                ],
                yerr=sn.lcs[control_index].t.loc[
                    good_ix, sn.lcs[control_index].colnames.dflux_new
                ],
                fmt="none",
                ecolor=SN_FLUX_COLORS[sn.filt],
                elinewidth=1,
                capsize=1.2,
                c=SN_FLUX_COLORS[sn.filt],
                alpha=0.5,
            )
            ax1.scatter(
                sn.lcs[control_index].t.loc[
                    good_ix, sn.lcs[control_index].colnames.mjd
                ],
                sn.lcs[control_index].t.loc[
                    good_ix, sn.lcs[control_index].colnames.flux
                ],
                s=marker_size,
                lw=marker_edgewidth,
                color=SN_FLUX_COLORS[sn.filt],
                marker="o",
                alpha=0.5,
                label="Cleaned measurements",
            )

            ax2.errorbar(
                sn.lcs[control_index].t.loc[
                    good_ix, sn.lcs[control_index].colnames.mjd
                ],
                sn.lcs[control_index].t.loc[
                    good_ix, sn.lcs[control_index].colnames.flux
                ],
                yerr=sn.lcs[control_index].t.loc[
                    good_ix, sn.lcs[control_index].colnames.dflux_new
                ],
                fmt="none",
                ecolor=SN_FLUX_COLORS[sn.filt],
                elinewidth=1,
                capsize=1.2,
                c=SN_FLUX_COLORS[sn.filt],
                alpha=0.5,
                zorder=5,
            )
            ax2.scatter(
                sn.lcs[control_index].t.loc[
                    good_ix, sn.lcs[control_index].colnames.mjd
                ],
                sn.lcs[control_index].t.loc[
                    good_ix, sn.lcs[control_index].colnames.flux
                ],
                s=marker_size,
                lw=marker_edgewidth,
                color=SN_FLUX_COLORS[sn.filt],
                marker="o",
                alpha=0.5,
                label="Cleaned measurements",
                zorder=5,
            )

        if sn.lcs[control_index].can_plot(bad_ix):
            ax2.errorbar(
                sn.lcs[control_index].t.loc[bad_ix, sn.lcs[control_index].colnames.mjd],
                sn.lcs[control_index].t.loc[
                    bad_ix, sn.lcs[control_index].colnames.flux
                ],
                yerr=sn.lcs[control_index].t.loc[
                    bad_ix, sn.lcs[control_index].colnames.dflux_new
                ],
                fmt="none",
                ecolor=SN_FLAGGED_FLUX_COLOR,
                elinewidth=1,
                capsize=1.2,
                c=SN_FLAGGED_FLUX_COLOR,
                alpha=0.5,
                zorder=10,
            )
            ax2.scatter(
                sn.lcs[control_index].t.loc[bad_ix, sn.lcs[control_index].colnames.mjd],
                sn.lcs[control_index].t.loc[
                    bad_ix, sn.lcs[control_index].colnames.flux
                ],
                s=marker_size,
                lw=marker_edgewidth,
                color=SN_FLAGGED_FLUX_COLOR,
                facecolors="none",
                edgecolors=SN_FLAGGED_FLUX_COLOR,
                marker="o",
                alpha=0.5,
                label="Flagged measurements",
                zorder=10,
            )

        lims = self.get_lims(sn=sn, custom_lims=custom_lims, flag=flag)
        if lims.get_xlims() is not None:
            ax1.set_xlim(lims.get_xlims())
            ax2.set_xlim(lims.get_xlims())
        if lims.get_ylims() is not None:
            ax1.set_ylim(lims.get_ylims())
            ax2.set_ylim(lims.get_ylims())

        ax1.legend(loc="upper right", facecolor="white", framealpha=1.0).set_zorder(100)
        ax2.legend(loc="upper right", facecolor="white", framealpha=1.0).set_zorder(100)

        if not save_filename is None:
            self.save_plot(save_filename)

        return fig

    def plot_cleaned_SN(
        self,
        sn: Supernova,
        flag: int,
        custom_lims: Optional[PlotLimits] = None,
        plot_controls: bool = True,
        plot_flagged: bool = False,
        save: bool = False,
        filename: str = "cleaned",
    ):
        fig, ax1 = plt.subplots(1, constrained_layout=True)
        fig.set_figwidth(7)
        fig.set_figheight(4)

        title = f"Cleaned SN {sn.tnsname}"
        if plot_controls and sn.num_controls > 0:
            title += f" & control light curves"
        title += f" {sn.filt}-band flux"
        ax1.set_title(title)

        ax1.minorticks_on()
        ax1.tick_params(direction="in", which="both")
        ax1.set_ylabel(r"Flux ($\mu$Jy)")
        ax1.set_xlabel(sn.colnames.mjd)
        ax1.axhline(linewidth=1, color="k")

        if plot_controls and sn.num_controls > 0:
            # plot control light curves
            label = f"Cleaned control measurements"
            for control_index in sn.get_control_lc_indices():
                lc = sn.lcs[control_index]
                good_ix = lc.get_good_indices(flag)

                if lc.can_plot(good_ix):
                    plt.errorbar(
                        lc.t.loc[good_ix, lc.colnames.mjd],
                        lc.t.loc[good_ix, lc.colnames.flux],
                        yerr=lc.t.loc[good_ix, lc.colnames.dflux_new],
                        fmt="none",
                        ecolor=CONTROL_FLUX_COLOR,
                        elinewidth=1.5,
                        capsize=1.2,
                        c=CONTROL_FLUX_COLOR,
                        alpha=0.5,
                        zorder=0,
                    )
                    plt.scatter(
                        lc.t.loc[good_ix, lc.colnames.mjd],
                        lc.t.loc[good_ix, lc.colnames.flux],
                        s=marker_size,
                        color=CONTROL_FLUX_COLOR,
                        marker="o",
                        alpha=0.5,
                        zorder=0,
                        label=label,
                    )

                if not label is None:
                    label = None

        sn_lc = sn.lcs[0]
        good_ix = sn_lc.get_good_indices(flag)

        if plot_flagged:
            bad_ix = sn_lc.get_bad_indices(flag)

            if sn_lc.can_plot(bad_ix):
                ax1.errorbar(
                    sn_lc.t.loc[bad_ix, sn_lc.colnames.mjd],
                    sn_lc.t.loc[bad_ix, sn_lc.colnames.flux],
                    yerr=sn_lc.t.loc[bad_ix, sn_lc.colnames.dflux_new],
                    fmt="none",
                    ecolor=SN_FLAGGED_FLUX_COLOR,
                    elinewidth=1,
                    capsize=1.2,
                    c=SN_FLAGGED_FLUX_COLOR,
                    alpha=0.5,
                    zorder=10,
                )
                ax1.scatter(
                    sn_lc.t.loc[bad_ix, sn_lc.colnames.mjd],
                    sn_lc.t.loc[bad_ix, sn_lc.colnames.flux],
                    s=marker_size,
                    lw=marker_edgewidth,
                    color=SN_FLAGGED_FLUX_COLOR,
                    facecolors="none",
                    edgecolors=SN_FLAGGED_FLUX_COLOR,
                    marker="o",
                    alpha=0.5,
                    label=f"Flagged SN {sn.tnsname} measurements",
                    zorder=10,
                )

        if sn_lc.can_plot(good_ix):
            plt.errorbar(
                sn_lc.t.loc[good_ix, sn_lc.colnames.mjd],
                sn_lc.t.loc[good_ix, sn_lc.colnames.flux],
                yerr=sn_lc.t.loc[good_ix, sn_lc.colnames.dflux_new],
                fmt="none",
                ecolor=SN_FLUX_COLORS[sn.filt],
                elinewidth=1,
                capsize=1.2,
                c=SN_FLUX_COLORS[sn.filt],
                alpha=0.5,
                zorder=10,
            )
            plt.scatter(
                sn_lc.t.loc[good_ix, sn_lc.colnames.mjd],
                sn_lc.t.loc[good_ix, sn_lc.colnames.flux],
                s=marker_size,
                lw=marker_edgewidth,
                color=SN_FLUX_COLORS[sn.filt],
                marker="o",
                alpha=0.5,
                zorder=10,
                label=f"Cleaned SN {sn.tnsname} measurements",
            )

        lims = self.get_lims(sn=sn, custom_lims=custom_lims, flag=flag)
        if lims.get_xlims() is not None:
            ax1.set_xlim(lims.get_xlims())
        if lims.get_ylims() is not None:
            ax1.set_ylim(lims.get_ylims())

        ax1.legend(loc="upper right", facecolor="white", framealpha=1.0).set_zorder(100)

        if save:
            self.save_plot(filename)

        return fig

    def plot_averaged_SN(
        self,
        avg_sn: AveragedSupernova,
        flag: int,
        custom_lims: Optional[PlotLimits] = None,
        plot_controls: bool = True,
        plot_flagged: bool = False,
        save: bool = False,
        filename: str = "averaged",
    ):
        fig, ax1 = plt.subplots(1, constrained_layout=True)
        fig.set_figwidth(7)
        fig.set_figheight(4)

        title = f"Cleaned & averaged SN {avg_sn.tnsname}"
        if plot_controls and avg_sn.num_controls > 0:
            title += f" & control light curves"
        title += f" {avg_sn.filt}-band flux"
        ax1.set_title(title)

        ax1.minorticks_on()
        ax1.tick_params(direction="in", which="both")
        ax1.set_ylabel(r"Flux ($\mu$Jy)")
        ax1.set_xlabel(avg_sn.colnames.mjd)
        ax1.axhline(linewidth=1, color="k")

        if plot_controls and avg_sn.num_controls > 0:
            # plot control light curves
            label = f"Cleaned & averaged control measurements"
            for control_index in avg_sn.get_control_lc_indices():
                lc = avg_sn.lcs[control_index]
                good_ix = lc.get_good_indices(flag)

                if lc.can_plot(good_ix):
                    plt.errorbar(
                        lc.t.loc[good_ix, lc.colnames.mjd],
                        lc.t.loc[good_ix, lc.colnames.flux],
                        yerr=lc.t.loc[good_ix, lc.colnames.dflux],
                        fmt="none",
                        ecolor=CONTROL_FLUX_COLOR,
                        elinewidth=1.5,
                        capsize=1.2,
                        c=CONTROL_FLUX_COLOR,
                        alpha=0.5,
                        zorder=0,
                    )
                    plt.scatter(
                        lc.t.loc[good_ix, lc.colnames.mjd],
                        lc.t.loc[good_ix, lc.colnames.flux],
                        s=marker_size,
                        color=CONTROL_FLUX_COLOR,
                        marker="o",
                        alpha=0.5,
                        zorder=0,
                        label=label,
                    )

                if not label is None:
                    label = None

        avg_sn_lc = avg_sn.lcs[0]
        good_ix = avg_sn_lc.get_good_indices(flag)

        if plot_flagged:
            bad_ix = avg_sn_lc.get_bad_indices(flag)

            if avg_sn_lc.can_plot(bad_ix):
                ax1.errorbar(
                    avg_sn_lc.t.loc[bad_ix, avg_sn_lc.colnames.mjd],
                    avg_sn_lc.t.loc[bad_ix, avg_sn_lc.colnames.flux],
                    yerr=avg_sn_lc.t.loc[bad_ix, avg_sn_lc.colnames.dflux],
                    fmt="none",
                    ecolor=SN_FLAGGED_FLUX_COLOR,
                    elinewidth=1,
                    capsize=1.2,
                    c=SN_FLAGGED_FLUX_COLOR,
                    alpha=0.5,
                    zorder=10,
                )
                ax1.scatter(
                    avg_sn_lc.t.loc[bad_ix, avg_sn_lc.colnames.mjd],
                    avg_sn_lc.t.loc[bad_ix, avg_sn_lc.colnames.flux],
                    s=marker_size,
                    lw=marker_edgewidth,
                    color=SN_FLAGGED_FLUX_COLOR,
                    facecolors="none",
                    edgecolors=SN_FLAGGED_FLUX_COLOR,
                    marker="o",
                    alpha=0.5,
                    label=f"Flagged averaged SN {avg_sn.tnsname} measurements",
                    zorder=10,
                )

        if avg_sn_lc.can_plot(good_ix):
            plt.errorbar(
                avg_sn_lc.t.loc[good_ix, avg_sn_lc.colnames.mjd],
                avg_sn_lc.t.loc[good_ix, avg_sn_lc.colnames.flux],
                yerr=avg_sn_lc.t.loc[good_ix, avg_sn_lc.colnames.dflux],
                fmt="none",
                ecolor=SN_FLUX_COLORS[avg_sn_lc.filt],
                elinewidth=1,
                capsize=1.2,
                c=SN_FLUX_COLORS[avg_sn_lc.filt],
                alpha=0.5,
                zorder=10,
            )
            plt.scatter(
                avg_sn_lc.t.loc[good_ix, avg_sn_lc.colnames.mjd],
                avg_sn_lc.t.loc[good_ix, avg_sn_lc.colnames.flux],
                s=marker_size,
                lw=marker_edgewidth,
                color=SN_FLUX_COLORS[avg_sn_lc.filt],
                marker="o",
                alpha=0.5,
                zorder=10,
                label=f"Cleaned & averaged SN {avg_sn.tnsname} measurements",
            )

        lims = self.get_lims(sn=avg_sn, custom_lims=custom_lims, flag=flag)
        if lims.get_xlims() is not None:
            ax1.set_xlim(lims.get_xlims())
        if lims.get_ylims() is not None:
            ax1.set_ylim(lims.get_ylims())

        ax1.legend(loc="upper right", facecolor="white", framealpha=1.0).set_zorder(100)

        if save:
            self.save_plot(filename)

        return fig

    def plot_limcuts(
        self,
        limcuts: LimCutsTable,
        cut: ChiSquareCut,
        save: bool = False,
        filename: str = "limcutstable",
    ):
        loss_color = "darkmagenta"
        contam_color = "teal"

        fig, ax1 = plt.subplots(1, constrained_layout=True)
        fig.set_figwidth(5.5)
        fig.set_figheight(3)

        ax1.set_title(f"Contamination and loss")

        ax1.minorticks_on()
        ax1.tick_params(direction="in", which="both")
        if cut.use_pre_mjd0_lc:
            ax1.set_ylabel(f"% pre-SN measurements")
        else:
            ax1.set_ylabel(f"% control measurements")
        ax1.set_xlabel("Chi-square cut")
        ax1.axhline(linewidth=1, color="k")

        ax1.plot(
            limcuts.t["PSF Chi-Square Cut"].values,
            limcuts.t["Ploss"].values,
            ms=3.5,
            color=loss_color,
            marker="o",
            label="Loss",
        )
        ax1.plot(
            limcuts.t["PSF Chi-Square Cut"].values,
            limcuts.t["Pcontamination"].values,
            ms=3.5,
            color=contam_color,
            marker="o",
            label="Contamination",
        )

        ax1.axvline(cut.max_value, color="k", linestyle="dashed", label="Selected cut")

        ax1.set_xlim(cut.min_cut, cut.max_cut)
        ax1.set_ylim(
            0, max(max(limcuts.t["Ploss"]), max(limcuts.t["Pcontamination"])) * 1.1
        )

        ax1.legend(
            facecolor="white", framealpha=1, bbox_to_anchor=(1.02, 1), loc="upper left"
        )

        if save:
            self.save_plot(filename)

        return fig

    def plot_uncert_est(
        self,
        sn: Supernova,
        custom_lims: PlotLimits,
        save: bool = False,
        filename: str = "uncert_est",
    ):
        lc = sn.lcs[0]
        if not f"{lc.colnames.dflux}_new" in lc.t.columns:
            print(
                f"WARNING: Cannot plot true uncertainties estimation due to missing {lc.colnames.dflux}_new column; skipping..."
            )
            return None

        fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True)
        fig.set_figwidth(7)
        fig.set_figheight(5)

        ax1.set_title(
            f"SN {sn.tnsname} {sn.filt}-band flux\nbefore true uncertainties estimation"
        )
        ax1.minorticks_on()
        ax1.tick_params(direction="in", which="both")
        ax1.get_xaxis().set_ticks([])
        ax1.set_ylabel(r"Flux ($\mu$Jy)")
        ax1.axhline(linewidth=1, color="k")

        ax2.set_title(f"after true uncertainties estimation")
        ax2.minorticks_on()
        ax2.tick_params(direction="in", which="both")
        ax2.set_ylabel(r"Flux ($\mu$Jy)")
        ax2.set_xlabel(lc.colnames.mjd)
        ax2.axhline(linewidth=1, color="k")

        ax1.errorbar(
            lc.t[lc.colnames.mjd],
            lc.t[lc.colnames.flux],
            yerr=lc.t[lc.colnames.dflux],
            fmt="none",
            ecolor=SN_FLUX_COLORS[lc.filt],
            elinewidth=1,
            capsize=1.2,
            c=SN_FLUX_COLORS[lc.filt],
            alpha=0.5,
        )
        ax1.scatter(
            lc.t[lc.colnames.mjd],
            lc.t[lc.colnames.flux],
            s=marker_size,
            lw=marker_edgewidth,
            color=SN_FLUX_COLORS[lc.filt],
            marker="o",
            alpha=0.5,
        )

        ax2.errorbar(
            lc.t[lc.colnames.mjd],
            lc.t[lc.colnames.flux],
            yerr=lc.t[lc.colnames.dflux_new],
            fmt="none",
            ecolor=SN_FLUX_COLORS[lc.filt],
            elinewidth=1,
            capsize=1.2,
            c=SN_FLUX_COLORS[lc.filt],
            alpha=0.5,
        )
        ax2.scatter(
            lc.t[lc.colnames.mjd],
            lc.t[lc.colnames.flux],
            s=marker_size,
            lw=marker_edgewidth,
            color=SN_FLUX_COLORS[lc.filt],
            marker="o",
            alpha=0.5,
        )

        lims = self.get_lims(sn=sn, custom_lims=custom_lims)
        if lims.get_xlims() is not None:
            ax1.set_xlim(lims.get_xlims())
            ax2.set_xlim(lims.get_xlims())
        if lims.get_ylims() is not None:
            ax1.set_ylim(lims.get_ylims())
            ax2.set_ylim(lims.get_ylims())

        if save:
            self.save_plot(filename)

        return fig

    def plot_template_correction(
        self,
        lc: LightCurve,
        custom_lims: PlotLimits,
        title=None,
        save: bool = False,
        filename: str = "template_correction",
    ):
        colors = ["salmon", "sandybrown", "darkseagreen"]

        region1_ix = lc.ix_inrange(lc.colnames.mjd, uplim=TEMPLATE_CHANGE_1_MJD)
        region2_ix = lc.ix_inrange(
            lc.colnames.mjd, lowlim=TEMPLATE_CHANGE_1_MJD, uplim=TEMPLATE_CHANGE_2_MJD
        )
        region3_ix = lc.ix_inrange(lc.colnames.mjd, lowlim=TEMPLATE_CHANGE_2_MJD)

        region1_mean = lc.get_mean(
            lc.colnames.flux, indices=region1_ix[-40:]
        )  # last 40 measurements before t1
        region2a_mean = lc.get_mean(
            lc.colnames.flux, indices=region2_ix[:40]
        )  # first 40 measurements after t1
        region2b_mean = lc.get_mean(
            lc.colnames.flux, indices=region2_ix[-40:]
        )  # last 40 measurements before t2
        region3_mean = lc.get_mean(
            lc.colnames.flux, indices=region3_ix[:40]
        )  # first 40 measurements after t2

        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], hspace=0.35, wspace=0.4)
        fig = plt.figure()
        fig.set_figwidth(6)
        fig.set_figheight(6)
        fig.tight_layout()

        # ax1: all template regions in different colors
        ax1 = plt.subplot(gs[0, :])
        if not title is None:
            ax1.set_title(title)
        ax1.axvline(
            x=TEMPLATE_CHANGE_1_MJD,
            color="k",
            linestyle="dotted",
            label="ATLAS template change",
            zorder=100,
        )
        ax1.axvline(x=TEMPLATE_CHANGE_2_MJD, color="k", linestyle="dotted", zorder=100)
        ax1.axhline(color="k", zorder=0)
        if custom_lims.get_xlims() is not None:
            ax1.set_xlim(custom_lims.get_xlims())
        if custom_lims.get_ylims() is not None:
            ax1.set_ylim(custom_lims.get_ylims())

        ax1.errorbar(
            lc.t.loc[region1_ix, lc.colnames.mjd],
            lc.t.loc[region1_ix, lc.colnames.flux],
            yerr=lc.t.loc[region1_ix, lc.colnames.dflux_new],
            fmt="none",
            ecolor=colors[0],
            elinewidth=1,
            capsize=1.2,
            c=SN_FLUX_COLORS[lc.filt],
            alpha=0.5,
            zorder=10,
        )
        ax1.scatter(
            lc.t.loc[region1_ix, lc.colnames.mjd],
            lc.t.loc[region1_ix, lc.colnames.flux],
            s=marker_size,
            lw=marker_edgewidth,
            color=colors[0],
            marker="o",
            alpha=0.5,
            zorder=10,
            label="Region 1 flux",
        )
        ax1.errorbar(
            lc.t.loc[region2_ix, lc.colnames.mjd],
            lc.t.loc[region2_ix, lc.colnames.flux],
            yerr=lc.t.loc[region2_ix, lc.colnames.dflux_new],
            fmt="none",
            ecolor=colors[1],
            elinewidth=1,
            capsize=1.2,
            c=SN_FLUX_COLORS[lc.filt],
            alpha=0.5,
            zorder=10,
        )
        ax1.scatter(
            lc.t.loc[region2_ix, lc.colnames.mjd],
            lc.t.loc[region2_ix, lc.colnames.flux],
            s=marker_size,
            lw=marker_edgewidth,
            color=colors[1],
            marker="o",
            alpha=0.5,
            zorder=10,
            label="Region 2 flux",
        )
        ax1.errorbar(
            lc.t.loc[region3_ix, lc.colnames.mjd],
            lc.t.loc[region3_ix, lc.colnames.flux],
            yerr=lc.t.loc[region3_ix, lc.colnames.dflux_new],
            fmt="none",
            ecolor=colors[2],
            elinewidth=1,
            capsize=1.2,
            c=SN_FLUX_COLORS[lc.filt],
            alpha=0.5,
            zorder=10,
        )
        ax1.scatter(
            lc.t.loc[region3_ix, lc.colnames.mjd],
            lc.t.loc[region3_ix, lc.colnames.flux],
            s=marker_size,
            lw=marker_edgewidth,
            color=colors[2],
            marker="o",
            alpha=0.5,
            zorder=10,
            label="Region 2 flux",
        )
        ax1.legend(
            facecolor="white", framealpha=1, loc="upper left", bbox_to_anchor=(1, 1)
        )

        # ax2: zoom in on first template change transition
        ax2 = plt.subplot(gs[1, 0])
        ax2.set_title("First template change", fontsize=12)
        ax2.axvline(x=TEMPLATE_CHANGE_1_MJD, color="k", linestyle="dotted", zorder=100)
        ax2.axhline(color="k", zorder=0)
        ax2.set_xlim(
            lc.t.loc[region1_ix[-40:][0], lc.colnames.mjd],
            lc.t.loc[region2_ix[:40][-1], lc.colnames.mjd],
        )
        if custom_lims.get_ylims() is not None:
            ax2.set_ylim(custom_lims.get_ylims())

        ax2.errorbar(
            lc.t.loc[region1_ix, lc.colnames.mjd],
            lc.t.loc[region1_ix, lc.colnames.flux],
            yerr=lc.t.loc[region1_ix, lc.colnames.dflux_new],
            fmt="none",
            ecolor=colors[0],
            elinewidth=1,
            capsize=1.2,
            c=SN_FLUX_COLORS[lc.filt],
            alpha=0.5,
            zorder=10,
        )
        ax2.scatter(
            lc.t.loc[region1_ix, lc.colnames.mjd],
            lc.t.loc[region1_ix, lc.colnames.flux],
            s=marker_size,
            lw=marker_edgewidth,
            color=colors[0],
            marker="o",
            alpha=0.5,
            zorder=10,
        )
        ax2.errorbar(
            lc.t.loc[region2_ix, lc.colnames.mjd],
            lc.t.loc[region2_ix, lc.colnames.flux],
            yerr=lc.t.loc[region2_ix, lc.colnames.dflux_new],
            fmt="none",
            ecolor=colors[1],
            elinewidth=1,
            capsize=1.2,
            c=SN_FLUX_COLORS[lc.filt],
            alpha=0.5,
            zorder=10,
        )
        ax2.scatter(
            lc.t.loc[region2_ix, lc.colnames.mjd],
            lc.t.loc[region2_ix, lc.colnames.flux],
            s=marker_size,
            lw=marker_edgewidth,
            color=colors[1],
            marker="o",
            alpha=0.5,
            zorder=10,
        )

        ax2.axhline(
            y=region1_mean, color=colors[0], linestyle="dashed", label="Region 1 mean"
        )
        ax2.axhline(
            y=region2a_mean, color=colors[1], linestyle="dashed", label="Region 2 mean"
        )
        ax2.legend(facecolor="white", framealpha=1)

        # ax3: zoom in on second template change transition
        ax3 = plt.subplot(gs[1, 1])
        ax3.set_title("Second template change", fontsize=12)
        ax3.axvline(x=TEMPLATE_CHANGE_2_MJD, color="k", linestyle="dotted", zorder=100)
        ax3.axhline(color="k", zorder=0)
        ax3.set_xlim(
            lc.t.loc[region2_ix[-40:][0], lc.colnames.mjd],
            lc.t.loc[region3_ix[:40][-1], lc.colnames.mjd],
        )
        if custom_lims.get_ylims() is not None:
            ax3.set_ylim(custom_lims.get_ylims())

        ax3.errorbar(
            lc.t.loc[region2_ix, lc.colnames.mjd],
            lc.t.loc[region2_ix, lc.colnames.flux],
            yerr=lc.t.loc[region2_ix, lc.colnames.dflux_new],
            fmt="none",
            ecolor=colors[1],
            elinewidth=1,
            capsize=1.2,
            c=SN_FLUX_COLORS[lc.filt],
            alpha=0.5,
            zorder=10,
        )
        ax3.scatter(
            lc.t.loc[region2_ix, lc.colnames.mjd],
            lc.t.loc[region2_ix, lc.colnames.flux],
            s=marker_size,
            lw=marker_edgewidth,
            color=colors[1],
            marker="o",
            alpha=0.5,
            zorder=10,
        )
        ax3.errorbar(
            lc.t.loc[region3_ix, lc.colnames.mjd],
            lc.t.loc[region3_ix, lc.colnames.flux],
            yerr=lc.t.loc[region3_ix, lc.colnames.dflux_new],
            fmt="none",
            ecolor=colors[2],
            elinewidth=1,
            capsize=1.2,
            c=SN_FLUX_COLORS[lc.filt],
            alpha=0.5,
            zorder=10,
        )
        ax3.scatter(
            lc.t.loc[region3_ix, lc.colnames.mjd],
            lc.t.loc[region3_ix, lc.colnames.flux],
            s=marker_size,
            lw=marker_edgewidth,
            color=colors[2],
            marker="o",
            alpha=0.5,
            zorder=10,
        )

        ax3.axhline(
            y=region2b_mean, color=colors[1], linestyle="dashed", label="Region 2 mean"
        )
        ax3.axhline(
            y=region3_mean, color=colors[2], linestyle="dashed", label="Region 3 mean"
        )
        ax3.legend(facecolor="white", framealpha=1)

        for ax in (ax1, ax2, ax3):
            ax.minorticks_on()
            ax.tick_params(direction="in", which="both")
            ax.set_xlabel("MJD")
            ax.set_ylabel(r"Flux ($\mu$Jy)")

        if save:
            self.save_plot(filename, bbox_inches="tight")

        return fig

    def plot_preSN(
        self,
        sn: Supernova,
        avg_sn: AveragedSupernova,
        flag: int,
        custom_lims: Optional[PlotLimits] = None,
        save: bool = False,
        filename: str = "pre_sn",
    ):
        fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True)
        fig.set_figwidth(4)
        fig.set_figheight(3.5)
        ax1.set_facecolor(FACE_COLOR)
        ax2.set_facecolor(FACE_COLOR)

        ax1.minorticks_on()
        ax1.set_xticklabels([])
        ax1.tick_params(direction="in", which="both")
        ax1.set_ylabel(r"Flux ($\mu$Jy)")
        ax1.axhline(linewidth=1.5, color="k", zorder=0)
        ax1.text(
            0.06,
            0.94,
            f"Pre-SN Light Curve",
            ha="left",
            va="top",
            transform=ax1.transAxes,
            fontsize=11,
            zorder=100,
        ).set_bbox(dict(facecolor=FACE_COLOR, alpha=0.8, edgecolor="silver"))

        ax2.minorticks_on()
        ax2.tick_params(direction="in", which="both")
        ax2.set_ylabel(r"Flux ($\mu$Jy)")
        ax2.set_xlabel("MJD")
        ax2.axhline(linewidth=1.5, color="k", zorder=0)
        ax2.text(
            0.06,
            0.94,
            f"Binned Pre-SN Light Curve",
            ha="left",
            va="top",
            transform=ax2.transAxes,
            fontsize=11,
            zorder=100,
        ).set_bbox(dict(facecolor=FACE_COLOR, alpha=0.8, edgecolor="silver"))

        # cleaned original light curve

        good_ix = sn.lcs[0].get_good_indices(flag)
        bad_ix = sn.lcs[0].get_bad_indices(flag)

        ax1.errorbar(
            sn.lcs[0].t.loc[good_ix, sn.colnames.mjd],
            sn.lcs[0].t.loc[good_ix, sn.colnames.flux],
            yerr=sn.lcs[0].t.loc[good_ix, sn.colnames.dflux_new],
            fmt="none",
            ecolor=SN_FLUX_COLORS[sn.filt],
            elinewidth=1.5,
            capsize=1.2,
            c=SN_FLUX_COLORS[sn.filt],
            alpha=0.5,
            zorder=0,
        )
        ax1.scatter(
            sn.lcs[0].t.loc[good_ix, sn.colnames.mjd],
            sn.lcs[0].t.loc[good_ix, sn.colnames.flux],
            s=marker_size,
            lw=marker_edgewidth,
            color=SN_FLUX_COLORS[sn.filt],
            marker="o",
            alpha=0.5,
            label=f"Cleaned Measurements",
            zorder=0,
        )

        ax1.errorbar(
            sn.lcs[0].t.loc[bad_ix, sn.colnames.mjd],
            sn.lcs[0].t.loc[bad_ix, sn.colnames.flux],
            yerr=sn.lcs[0].t.loc[bad_ix, sn.colnames.dflux_new],
            fmt="none",
            ecolor=SN_FLAGGED_FLUX_COLOR,
            elinewidth=1.5,
            capsize=1.2,
            c=SN_FLAGGED_FLUX_COLOR,
            alpha=0.5,
            zorder=10,
        )
        ax1.scatter(
            sn.lcs[0].t.loc[bad_ix, sn.colnames.mjd],
            sn.lcs[0].t.loc[bad_ix, sn.colnames.flux],
            s=marker_size,
            lw=marker_edgewidth,
            facecolors="none",
            edgecolors=SN_FLAGGED_FLUX_COLOR,
            marker="o",
            alpha=0.5,
            label=f"Flagged Measurements",
            zorder=10,
        )

        # averaged light curve

        good_ix = avg_sn.lcs[0].get_good_indices(flag)
        bad_ix = avg_sn.lcs[0].get_bad_indices(flag)

        ax2.errorbar(
            avg_sn.lcs[0].t.loc[good_ix, avg_sn.colnames.mjd],
            avg_sn.lcs[0].t.loc[good_ix, avg_sn.colnames.flux],
            yerr=avg_sn.lcs[0].t.loc[good_ix, avg_sn.colnames.dflux],
            fmt="none",
            ecolor=SN_FLUX_COLORS[sn.filt],
            elinewidth=1.5,
            capsize=1.2,
            c=SN_FLUX_COLORS[sn.filt],
            alpha=0.5,
            zorder=0,
        )
        ax2.scatter(
            avg_sn.lcs[0].t.loc[good_ix, avg_sn.colnames.mjd],
            avg_sn.lcs[0].t.loc[good_ix, avg_sn.colnames.flux],
            s=marker_size,
            lw=marker_edgewidth,
            color=SN_FLUX_COLORS[sn.filt],
            marker="o",
            alpha=0.5,
            label=f"Cleaned Measurements",
            zorder=0,
        )

        ax2.errorbar(
            avg_sn.lcs[0].t.loc[bad_ix, avg_sn.colnames.mjd],
            avg_sn.lcs[0].t.loc[bad_ix, avg_sn.colnames.flux],
            yerr=avg_sn.lcs[0].t.loc[bad_ix, avg_sn.colnames.dflux],
            fmt="none",
            ecolor=SN_FLAGGED_FLUX_COLOR,
            elinewidth=1.5,
            capsize=1.2,
            c=SN_FLAGGED_FLUX_COLOR,
            alpha=0.5,
            zorder=10,
        )
        ax2.scatter(
            avg_sn.lcs[0].t.loc[bad_ix, avg_sn.colnames.mjd],
            avg_sn.lcs[0].t.loc[bad_ix, avg_sn.colnames.flux],
            s=marker_size,
            lw=marker_edgewidth,
            facecolors="none",
            edgecolors=SN_FLAGGED_FLUX_COLOR,
            marker="o",
            alpha=0.5,
            label=f"Flagged Measurements",
            zorder=10,
        )

        ax2.legend(
            facecolor="white",
            edgecolor="silver",
            fontsize=9,
            framealpha=0.8,
            handletextpad=0.1,
            loc="lower right",
            borderaxespad=1,
            ncol=1,
        ).set_zorder(100)

        lims = self.get_lims(sn=sn, custom_lims=custom_lims, flag=flag)
        if lims.get_xlims() is not None:
            ax1.set_xlim(lims.get_xlims())
            ax2.set_xlim(lims.get_xlims())
        if lims.get_ylims() is not None:
            ax1.set_ylim(lims.get_ylims())
            ax2.set_ylim(lims.get_ylims())

        if save:
            self.save_plot(filename, bbox_inches="tight")

        """
        if plot_mjd_ranges and not mjd_ranges is None:
            # valid mjd ranges
            for mjd_range in mjd_ranges:
                ax1.axvline(mjd_range[0], color="k", linestyle="dashed", zorder=100)
                ax1.axvline(mjd_range[1], color="k", linestyle="dashed", zorder=100)
                ax2.axvline(mjd_range[0], color="k", linestyle="dashed", zorder=100)
                ax2.axvline(mjd_range[1], color="k", linestyle="dashed", zorder=100)
        """


class PlotPdf(Plot):
    def __init__(self, output_dir, tnsname, filt="o"):
        Plot.__init__(self)
        self.filename = f"{output_dir}/{tnsname}.{filt}.plots.pdf"
        self.pdf = PdfPages(self.filename)

    def save_pdf(self):
        print("\nSaving PDF of plots...\n")
        self.pdf.close()

    def plot_SN(
        self,
        sn: Supernova,
        custom_lims: Optional[PlotLimits] = None,
        plot_controls: bool = True,
        plot_template_changes: bool = True,
        save: bool = False,
        filename: str = "original",
    ):
        print(
            f'Plotting original SN{" and control light curves" if plot_controls else ""}...'
        )
        fig = super().plot_SN(
            sn, custom_lims, plot_controls, plot_template_changes, save, filename
        )
        self.pdf.savefig(fig)

    def plot_cut(
        self,
        sn: Supernova,
        flag: int,
        custom_lims: Optional[PlotLimits] = None,
        title: str | None = None,
        save_filename: str = None,
    ):
        print(f"Plotting cut for flag {hex(flag)}...")
        fig = super().plot_cut(sn, flag, custom_lims, title, save_filename)
        self.pdf.savefig(fig)

    def plot_cleaned_SN(
        self,
        sn: Supernova,
        flag: int,
        custom_lims: Optional[PlotLimits] = None,
        plot_controls: bool = True,
        plot_flagged: bool = True,
        save: bool = False,
        filename: str = "cleaned",
    ):
        print(
            f'Plotting cleaned SN{" and control light curves" if plot_controls else ""} using flag {hex(flag)}...'
        )
        fig = super().plot_cleaned_SN(
            sn, flag, custom_lims, plot_controls, plot_flagged, save, filename
        )
        self.pdf.savefig(fig)

    def plot_averaged_SN(
        self,
        avg_sn: AveragedSupernova,
        flag: int,
        custom_lims: Optional[PlotLimits] = None,
        plot_controls: bool = True,
        plot_flagged: bool = True,
        save: bool = False,
        filename: str = "averaged",
    ):
        print(
            f'Plotting averaged SN{" and control light curves" if plot_controls else ""} using flag {hex(flag)}...'
        )
        fig = super().plot_averaged_SN(
            avg_sn, flag, custom_lims, plot_controls, plot_flagged, save, filename
        )
        self.pdf.savefig(fig)

    def plot_limcuts(
        self,
        limcuts: LimCutsTable,
        cut: ChiSquareCut,
        save: bool = False,
        filename: str = "limcutstable",
    ):
        print("Plotting LimCutsTable...")
        fig = super().plot_limcuts(limcuts, cut, save, filename)
        self.pdf.savefig(fig)

    def plot_uncert_est(
        self,
        sn: Supernova,
        custom_lims: Optional[PlotLimits] = None,
        save: bool = False,
        filename: str = "uncert_est",
    ):
        print("Plotting true uncertainties estimation...")
        fig = super().plot_uncert_est(sn, custom_lims, save, filename)
        if not fig is None:
            self.pdf.savefig(fig)

    def plot_template_correction(self, lc: LightCurve):
        print("Plotting ATLAS template chanages correction...")
        fig = super().plot_template_correction(lc)
        self.pdf.savefig(fig)
