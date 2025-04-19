#!/usr/bin/env python

import os
from typing import Dict, List, Optional
import matplotlib
from matplotlib import gridspec
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from lightcurve import (
    LimCutsTable,
    LightCurve,
    SimDetecLightCurve,
    SimDetecSupernova,
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
MARKER_SIZE = 30
MARKER_EDGEWIDTH = 1.5

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
        sn: Supernova,
        control_index: int = 0,
        flag: Optional[int] = None,
        custom_lims: Optional[PlotLimits] = None,
        indices: Optional[List[int]] = None,
        pre_sn: bool = False,
    ) -> PlotLimits:
        lims = PlotLimits()
        mjd0 = sn.mjd0 if pre_sn else None

        # get auto xlims
        lims.set_xlims(sn.lcs[control_index].get_xlims(mjd0=mjd0))
        # override auto xlims with custom ones where they exist
        if custom_lims is not None and custom_lims.get_xlims() is not None:
            lims.set_xlims(custom_lims.get_xlims())

        # get auto ylims
        lims.set_ylims(
            sn.lcs[control_index].get_ylims(indices=indices, flag=flag, mjd0=mjd0)
        )
        # override auto ylims with custom ones where they exist
        if custom_lims is not None and custom_lims.get_ylims() is not None:
            lims.set_ylims(custom_lims.get_ylims())

        return lims

    def get_snr_lims(
        self,
        sn: SimDetecSupernova,
        all_fom: pd.Series,
        fom_limit: Optional[float] = None,
        custom_lims: Optional[PlotLimits] = None,
    ) -> PlotLimits:
        lims = PlotLimits()

        # get auto xlims using min and max of mjd ranges
        lims.set_xlims([sn._mjd_ranges[0][0], sn._mjd_ranges[-1][-1]])

        # override auto xlims with custom ones where they exist
        if custom_lims is not None and custom_lims.get_xlims() is not None:
            lims.set_xlims(custom_lims.get_xlims())

        # get auto ylims using min and max of FOM distribution
        lims.set_ylims(
            [
                min(all_fom) * 1.5,
                (
                    max(max(all_fom) * 1.2, fom_limit * 1.4)
                    if fom_limit
                    else max(all_fom) * 1.2
                ),
            ]
        )

        # override auto ylims with custom ones where they exist
        if custom_lims is not None and custom_lims.get_ylims() is not None:
            lims.set_ylims(custom_lims.get_ylims())

        return lims

    def _setup_ax(
        self,
        ax: Axes,
        lims: PlotLimits,
        xlabel: bool = True,
        ylabel: bool = True,
        xticks: bool = True,
        yticks: bool = True,
        axhline: bool = True,
    ):
        ax.minorticks_on()
        ax.tick_params(direction="in", which="both")
        ax.set_facecolor(FACE_COLOR)
        if axhline:
            ax.axhline(color="k", linewidth=1.5, zorder=0)

        if lims is not None and lims.get_xlims() is not None:
            ax.set_xlim(lims.get_xlims())
        if lims is not None and lims.get_ylims() is not None:
            ax.set_ylim(lims.get_ylims())

        if xlabel:
            ax.set_xlabel("MJD")
        if not xticks:
            ax.set_xticklabels([])

        if ylabel:
            ax.set_ylabel(r"Flux ($\mu$Jy)")
        if not yticks:
            ax.set_yticklabels([])

    def _plot_snr(
        self,
        ax: Axes,
        obj: SimDetecSupernova | SimDetecLightCurve,
        control_index: int,
        color: str,
        y_colname_attr: Optional[str] = "snrsumnorm",
        indices: Optional[List[int]] = None,
        label: Optional[str] = None,
    ):
        if isinstance(obj, SimDetecSupernova):
            obj = obj.lcs[control_index]

        if indices is None:
            indices = obj.getindices()

        y_colname = getattr(obj.colnames, y_colname_attr)
        if not obj.can_plot(indices, columns=[y_colname]):
            print(
                f"WARNING: Light curve (control index #{control_index}) '{y_colname_attr}' column cannot be plotted with indices of length {len(indices)}; skipping..."
            )
            return

        ax.plot(
            obj.t.loc[indices, obj.colnames.mjdbin],
            obj.t.loc[indices, y_colname],
            color=color,
            linewidth=1.5,
            alpha=0.7,
            label=label,
        )

    def _plot_lc(
        self,
        ax: Axes,
        obj: Supernova | LightCurve,
        control_index: int,
        color: str,
        indices: Optional[List[int]] = None,
        label: Optional[str] = None,
        open: bool = False,
    ):
        if isinstance(obj, Supernova):
            obj = obj.lcs[control_index]

        if indices is None:
            indices = obj.getindices()
        if not obj.can_plot(indices):
            print(
                f"WARNING: Light curve (control index #{control_index}) cannot be plotted with indices of length {len(indices)}; skipping..."
            )
            return

        ax.errorbar(
            obj.t.loc[indices, obj.colnames.mjd],
            obj.t.loc[indices, obj.colnames.flux],
            yerr=obj.t.loc[indices, obj.colnames.dflux_new],
            fmt="none",
            ecolor=color,
            elinewidth=1.5,
            capsize=1.2,
            c=color,
            alpha=0.5,
            zorder=10,
        )
        ax.scatter(
            obj.t.loc[indices, obj.colnames.mjd],
            obj.t.loc[indices, obj.colnames.flux],
            s=MARKER_SIZE,
            lw=MARKER_EDGEWIDTH,
            color=color,
            marker="o",
            alpha=0.5,
            zorder=10,
            label=label,
            facecolors="none" if open else None,
            edgecolors=color if open else None,
        )

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

        lims = self.get_lims(sn, custom_lims=custom_lims)
        self._setup_ax(ax1, lims)

        title = f"SN {sn.tnsname}"
        if plot_controls and sn.num_controls > 0:
            title += f" & control light curves"
        title += f" {sn.filt}-band flux"
        ax1.set_title(title)

        if plot_controls and sn.num_controls > 0:
            # plot control light curves
            label = f"{sn.num_controls} control light curves"
            for control_index in sn.control_lc_indices:
                self._plot_lc(ax1, sn, control_index, CONTROL_FLUX_COLOR, label=label)
                if not label is None:
                    label = None

        preMJD0_ix = sn.lcs[0].get_preMJD0_indices(sn.mjd0)
        postMJD0_ix = sn.lcs[0].get_postMJD0_indices(sn.mjd0)

        # plot pre-MJD0 SN light curve
        self._plot_lc(
            ax1, sn, 0, "magenta", indices=preMJD0_ix, label="Pre-MJD0 light curve"
        )

        # plot post-MJD0 SN light curve
        self._plot_lc(
            ax1, sn, 0, "lime", indices=postMJD0_ix, label="Post-MJD0 light curve"
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
        fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True)
        ax1: Axes
        ax2: Axes
        fig.set_figwidth(7)
        fig.set_figheight(5)

        if not title:
            title = "Cut"
        fig.suptitle(f"{title} (flag {hex(flag)})")

        lims = self.get_lims(sn, custom_lims=custom_lims, flag=flag)
        self._setup_ax(ax1, lims, ylabel=False, xticks=False, xlabel=False)
        self._setup_ax(ax2, lims, ylabel=False)

        good_ix = sn.lcs[control_index].get_good_indices(flag)
        bad_ix = sn.lcs[control_index].get_bad_indices(flag)

        self._plot_lc(
            ax1,
            sn,
            control_index,
            SN_FLUX_COLORS[sn.filt],
            indices=good_ix,
            label="Cleaned measurements",
        )
        self._plot_lc(
            ax2,
            sn,
            control_index,
            SN_FLUX_COLORS[sn.filt],
            indices=good_ix,
            label="Cleaned measurements",
        )

        self._plot_lc(
            ax1,
            sn,
            control_index,
            SN_FLAGGED_FLUX_COLOR,
            indices=bad_ix,
            label="Flagged measurements",
            open=True,
        )

        fig.supylabel(r"Flux ($\mu$Jy)")

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
        ax1: Axes
        fig.set_figwidth(7)
        fig.set_figheight(4)

        title = f"Cleaned SN {sn.tnsname}"
        if plot_controls and sn.num_controls > 0:
            title += f" & control light curves"
        title += f" {sn.filt}-band flux"
        ax1.set_title(title)

        lims = self.get_lims(sn, custom_lims=custom_lims, flag=flag)
        self._setup_ax(ax1, lims)

        if plot_controls and sn.num_controls > 0:
            # plot control light curves
            label = f"Cleaned control measurements"
            for control_index in sn.control_lc_indices:
                good_ix = sn.lcs[control_index].get_good_indices(flag)
                self._plot_lc(
                    ax1,
                    sn,
                    control_index,
                    CONTROL_FLUX_COLOR,
                    indices=good_ix,
                    label=label,
                )
                if not label is None:
                    label = None

        if plot_flagged:
            bad_ix = sn.lcs[0].get_bad_indices(flag)
            self._plot_lc(
                ax1,
                sn,
                0,
                SN_FLAGGED_FLUX_COLOR,
                indices=bad_ix,
                label=f"Flagged SN measurements",
                open=True,
            )

        good_ix = sn.lcs[0].get_good_indices(flag)
        self._plot_lc(
            ax1,
            sn,
            0,
            SN_FLUX_COLORS[sn.filt],
            indices=good_ix,
            label=f"Cleaned SN measurements",
        )

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
        ax1: Axes
        fig.set_figwidth(7)
        fig.set_figheight(4)

        title = f"Cleaned & averaged SN {avg_sn.tnsname}"
        if plot_controls and avg_sn.num_controls > 0:
            title += f" & control light curves"
        title += f" {avg_sn.filt}-band flux"
        ax1.set_title(title)

        lims = self.get_lims(avg_sn, custom_lims=custom_lims, flag=flag)
        self._setup_ax(ax1, lims)

        if plot_controls and avg_sn.num_controls > 0:
            # plot control light curves
            label = f"Cleaned control bins"
            for control_index in avg_sn.control_lc_indices:
                good_ix = avg_sn.lcs[control_index].get_good_indices(flag)
                self._plot_lc(
                    ax1,
                    avg_sn,
                    control_index,
                    CONTROL_FLUX_COLOR,
                    indices=good_ix,
                    label=label,
                )

                if not label is None:
                    label = None

        if plot_flagged:
            bad_ix = avg_sn.lcs[0].get_bad_indices(flag)
            self._plot_lc(
                ax1,
                avg_sn,
                0,
                SN_FLAGGED_FLUX_COLOR,
                indices=bad_ix,
                label=f"Flagged SN bins",
                open=True,
            )

        good_ix = avg_sn.lcs[control_index].get_good_indices(flag)
        self._plot_lc(
            ax1,
            avg_sn,
            0,
            SN_FLUX_COLORS[avg_sn.filt],
            indices=good_ix,
            label=f"Cleaned SN bins",
        )

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
        ax1: Axes
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
        ax1: Axes
        ax2: Axes
        fig.set_figwidth(7)
        fig.set_figheight(5)

        lims = self.get_lims(sn, custom_lims=custom_lims)

        ax1.set_title(
            f"SN {sn.tnsname} {sn.filt}-band flux\nbefore true uncertainties estimation"
        )
        self._setup_ax(ax1, lims, xticks=False, xlabel=False)

        ax2.set_title(f"after true uncertainties estimation")
        self._setup_ax(ax2, lims)

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
            s=MARKER_SIZE,
            lw=MARKER_EDGEWIDTH,
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
            s=MARKER_SIZE,
            lw=MARKER_EDGEWIDTH,
            color=SN_FLUX_COLORS[lc.filt],
            marker="o",
            alpha=0.5,
        )

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

        # last 40 measurements before t1
        region1_mean = lc.get_mean(lc.colnames.flux, indices=region1_ix[-40:])
        # first 40 measurements after t1
        region2a_mean = lc.get_mean(lc.colnames.flux, indices=region2_ix[:40])
        # last 40 measurements before t2
        region2b_mean = lc.get_mean(lc.colnames.flux, indices=region2_ix[-40:])
        # first 40 measurements after t2
        region3_mean = lc.get_mean(lc.colnames.flux, indices=region3_ix[:40])

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

        self._plot_lc(ax1, lc, 0, colors[0], region1_ix, label="Region 1 flux")
        self._plot_lc(ax1, lc, 0, colors[1], region2_ix, label="Region 2 flux")
        self._plot_lc(ax1, lc, 0, colors[2], region3_ix, label="Region 3 flux")

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

        self._plot_lc(ax2, lc, 0, colors[0], region1_ix)
        self._plot_lc(ax2, lc, 0, colors[1], region2_ix)

        ax2.axhline(
            y=region1_mean, color=colors[0], linestyle="dashed", label="Region 1 mean"
        )
        ax2.axhline(
            y=region2a_mean, color=colors[1], linestyle="dashed", label="Region 2 mean"
        )

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

        self._plot_lc(ax3, lc, 0, colors[1], region2_ix)
        self._plot_lc(ax3, lc, 0, colors[2], region3_ix)

        ax3.axhline(
            y=region2b_mean, color=colors[1], linestyle="dashed", label="Region 2 mean"
        )
        ax3.axhline(
            y=region3_mean, color=colors[2], linestyle="dashed", label="Region 3 mean"
        )

        ax1.legend(
            facecolor="white", framealpha=1, loc="upper left", bbox_to_anchor=(1, 1)
        )
        ax2.legend(facecolor="white", framealpha=1)
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
        ax1: Axes
        ax2: Axes
        fig.set_figwidth(4)
        fig.set_figheight(3.5)

        lims = self.get_lims(sn, custom_lims=custom_lims, flag=flag, pre_sn=True)

        self._setup_ax(ax1, lims, xlabel=False, xticks=False)
        ax1.set_title("Pre-SN Light Curve", fontsize=12)

        self._setup_ax(ax2, lims)
        ax2.set_title("Binned Pre-SN Light Curve", fontsize=12)

        # cleaned original light curve
        self._plot_lc(
            ax1,
            sn,
            0,
            SN_FLUX_COLORS[sn.filt],
            indices=sn.lcs[0].get_good_indices(flag),
            label="Cleaned",
        )
        self._plot_lc(
            ax1,
            sn,
            0,
            SN_FLAGGED_FLUX_COLOR,
            indices=sn.lcs[0].get_bad_indices(flag),
            label="Flagged",
            open=True,
        )

        # averaged light curve
        self._plot_lc(
            ax2,
            avg_sn,
            0,
            SN_FLUX_COLORS[sn.filt],
            indices=avg_sn.lcs[0].get_good_indices(flag),
            label="Cleaned",
        )
        self._plot_lc(
            ax2,
            avg_sn,
            0,
            SN_FLAGGED_FLUX_COLOR,
            indices=avg_sn.lcs[0].get_bad_indices(flag),
            label="Flagged",
            open=True,
        )

        ax2.legend(
            facecolor="white",
            edgecolor="silver",
            fontsize=9,
            framealpha=0.8,
            handletextpad=0.1,
            loc="upper left",
            borderaxespad=1,
            ncol=1,
        ).set_zorder(100)

        if save:
            self.save_plot(filename, bbox_inches="tight")

        return fig

    def plot_mjd_ranges(
        self,
        sn: Supernova,
        avg_sn: AveragedSupernova,
        flag: int,
        mjd_ranges: List[List],
        custom_lims: Optional[PlotLimits] = None,
        suptitle: Optional[str] = None,
        range_color: str = "gray",
        save: bool = False,
        filename: str = "mjd_ranges",
    ):

        fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True)
        ax1: Axes
        ax2: Axes
        fig.set_figwidth(4)
        fig.set_figheight(3.5)
        if suptitle:
            fig.suptitle(suptitle)

        lims = self.get_lims(sn, custom_lims=custom_lims, flag=flag)

        self._setup_ax(ax1, lims, xlabel=False, xticks=False)
        ax1.set_title("SN Light Curve", fontsize=12)

        self._setup_ax(ax2, lims)
        ax2.set_title("Binned SN Light Curve", fontsize=12)

        # cleaned original light curve
        self._plot_lc(
            ax1,
            sn,
            0,
            SN_FLUX_COLORS[sn.filt],
            indices=sn.lcs[0].get_good_indices(flag),
            label="Cleaned",
        )
        self._plot_lc(
            ax1,
            sn,
            0,
            SN_FLAGGED_FLUX_COLOR,
            indices=sn.lcs[0].get_bad_indices(flag),
            label="Flagged",
            open=True,
        )

        # averaged light curve
        self._plot_lc(
            ax2,
            avg_sn,
            0,
            SN_FLUX_COLORS[sn.filt],
            indices=avg_sn.lcs[0].get_good_indices(flag),
            label="Cleaned",
        )
        self._plot_lc(
            ax2,
            avg_sn,
            0,
            SN_FLAGGED_FLUX_COLOR,
            indices=avg_sn.lcs[0].get_bad_indices(flag),
            label="Flagged",
            open=True,
        )

        if mjd_ranges is not None:
            for ax in [ax1, ax2]:
                for mjd_range in mjd_ranges:
                    ax.axvspan(
                        mjd_range[0],
                        mjd_range[1],
                        color=range_color,
                        alpha=0.2,
                        zorder=0,
                    )

        ax2.legend(
            facecolor="white",
            edgecolor="silver",
            fontsize=9,
            framealpha=0.8,
            handletextpad=0.1,
            loc="upper left",
            borderaxespad=1,
            ncol=1,
        ).set_zorder(100)

        if save:
            self.save_plot(filename, bbox_inches="tight")

        return fig

    def plot_all_controls(
        self,
        sn: Supernova,
        flag: int,
        custom_lims: Optional[PlotLimits] = None,
        two_columns: bool = False,
        save: bool = False,
        filename: str = "all_controls",
    ):
        if sn.num_controls < 1:
            print("WARNING: No control light curves to plot")
            return

        if two_columns and sn.num_controls % 2 != 0:
            raise RuntimeError(
                f"Number of control light curves ({sn.num_controls}) must be even; "
                "set two_columns=False for one column"
            )

        control_indices = sn.control_lc_indices
        lims = self.get_lims(
            sn, control_index=control_indices[0], custom_lims=custom_lims
        )

        # set up figure and axes
        if two_columns:
            num_rows = sn.num_controls // 2
            fig, axes = plt.subplots(
                num_rows, 2, constrained_layout=True, figsize=(7, num_rows * 1.2)
            )
            axes = axes.flatten()
        else:
            fig, axes = plt.subplots(
                sn.num_controls,
                1,
                constrained_layout=True,
                figsize=(5, sn.num_controls),
            )
            axes = np.atleast_1d(axes)

        # loop over control light curves
        for idx, control_index in enumerate(control_indices):
            ax: Axes = axes[idx]

            is_last_row = idx >= len(control_indices) - (2 if two_columns else 1)
            is_rightmost_col = idx % 2 == 1 if two_columns else False
            self._setup_ax(
                ax,
                lims,
                xlabel=is_last_row,
                xticks=is_last_row,
                ylabel=False,
                yticks=not is_rightmost_col,
            )

            good_ix = sn.lcs[control_index].get_good_indices(flag)
            self._plot_lc(
                ax, sn, control_index, SELECT_CONTROL_FLUX_COLOR, indices=good_ix
            )

            label_text = (
                f"Binned & Cleaned Control Light Curve #{control_index}"
                if idx == 0
                else f"#{control_index}"
            )
            ax.text(
                0.03,
                0.92,
                label_text,
                ha="left",
                va="top",
                transform=ax.transAxes,
                fontsize=11,
                zorder=20,
            )

        fig.supylabel(r"Flux (µJy)")

        if save:
            self.save_plot(filename, bbox_inches="tight")

        return fig

    def plot_binned_examples(
        self,
        avg_sn: AveragedSupernova,
        select_control_index: int,
        flag: int,
        custom_lims: Optional[PlotLimits] = None,
        save: bool = False,
        filename: str = "binned_examples",
    ):
        fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True)
        ax1: Axes
        ax2: Axes
        fig.set_figwidth(4)
        fig.set_figheight(3.5)

        lims = self.get_lims(avg_sn, custom_lims=custom_lims, flag=flag, pre_sn=True)

        self._setup_ax(ax1, lims, xlabel=False, xticks=False)
        ax1.set_title("Binned & Cleaned Pre-SN Light Curve", fontsize=12)
        self._plot_lc(
            ax1,
            avg_sn,
            0,
            SN_FLUX_COLORS[avg_sn.filt],
            indices=avg_sn.lcs[0].get_good_indices(flag),
        )

        self._setup_ax(ax2, lims)
        ax2.set_title(
            f"Binned & Cleaned Control Light Curve #{select_control_index}", fontsize=12
        )
        self._plot_lc(
            ax2,
            avg_sn,
            select_control_index,
            SELECT_CONTROL_FLUX_COLOR,
            indices=avg_sn.lcs[select_control_index].get_good_indices(flag),
        )

        if save:
            self.save_plot(filename, bbox_inches="tight")

        return fig

    def plot_fom(
        self,
        sn: SimDetecSupernova,
        all_fom_dict: Dict[int, pd.Series],
        sigma_kerns: List[int],
        select_control_index: int,
        flag: int,
        fom_limits: Optional[Dict[int, float]] = None,
        save: bool = False,
        filename: str = "all_fom",
    ):
        fig, axes = plt.subplots(
            len(sigma_kerns), 2, gridspec_kw={"width_ratios": [4, 1]}
        )
        fig.set_figheight(len(sigma_kerns) * 1.9)
        fig.set_figwidth(5.5)
        fig.subplots_adjust(wspace=0.07, hspace=0.07, bottom=0.1, top=0.9)

        max_freq = 0
        for sigma_kern, all_fom in all_fom_dict.items():
            counts, _ = np.histogram(all_fom.dropna(), bins=20)
            max_freq = max(max_freq, counts.max())

        for i, row in enumerate(axes):
            sn.apply_rolling_sums(sigma_kern, flag=flag, pre_mjd0_ix=True)

            if len(sigma_kerns) == 1:
                sigma_kern = sigma_kerns[0]
                ax1, ax2 = axes[0], axes[1]
            else:
                sigma_kern = sigma_kerns[i]
                ax1, ax2 = axes[i][0], axes[i][1]
            ax1: Axes
            ax2: Axes

            lims = self.get_snr_lims(
                sn,
                all_fom_dict[sigma_kern],
                fom_limit=fom_limits[sigma_kern] if fom_limits else None,
            )
            self._setup_ax(ax1, lims, ylabel=False)
            self._setup_ax(ax2, lims, ylabel=False, yticks=False, axhline=False)
            ax1.set_ylabel(r"$\Sigma_{\rm FOM}$")
            if i >= len(sigma_kerns) - 1:  # bottom row
                ax1.set_xlabel("MJD")
                ax2.set_xlabel("Freq")
            else:
                ax1.set_xticklabels([])
                ax2.set_xticklabels([])

            # control lc fom
            label_control_lc_indices = [
                x for x in sn.control_lc_indices if x != select_control_index
            ]
            for control_index in sn.control_lc_indices:
                if control_index == select_control_index:
                    continue
                label = None
                if control_index == sn.control_lc_indices[0]:
                    label = f"{len(sn.control_lc_indices) - 1} Control Light Curves (#s: {label_control_lc_indices})"
                self._plot_snr(ax1, sn, control_index, CONTROL_FOM_COLOR, label=label)

            # selected control lc fom
            self._plot_snr(
                ax1,
                sn,
                select_control_index,
                SELECT_CONTROL_FOM_COLOR,
                label=f"Selected Control Light Curve #{select_control_index}",
            )

            # pre-SN lc fom
            self._plot_snr(ax1, sn, 0, SN_FOM_COLOR, label="Pre-SN Light Curve")

            # sigma_kern label
            ax1.text(
                0.98,
                0.07,
                r"$\sigma_{\rm kernel}$ = " + str(sigma_kern),
                ha="right",
                va="bottom",
                transform=ax1.transAxes,
                fontsize=11,
                zorder=40,
            )

            # fom limit
            if fom_limits:
                ax1.axhline(
                    fom_limits[sigma_kern],
                    linewidth=1.5,
                    color="k",
                    linestyle=FOM_LIMIT_LS,
                    zorder=40,
                )
                ax2.axhline(
                    fom_limits[sigma_kern],
                    linewidth=1.5,
                    color="k",
                    linestyle=FOM_LIMIT_LS,
                    zorder=40,
                )
                ax1.text(
                    0.05,
                    1.1 * fom_limits[sigma_kern],
                    r"$\Sigma_{\rm FOM, limit}$ = " + str(fom_limits[sigma_kern]),
                    color="k",
                    transform=ax1.get_yaxis_transform(),
                    zorder=40,
                )

            # fom distribution
            ax2.set_ylim(ax1.get_ylim())
            ax2.set_xlim(0, max_freq * 1.1)
            ax2.hist(
                all_fom_dict[sigma_kern],
                bins=20,
                orientation="horizontal",
                color=CONTROL_FOM_COLOR,
            )

            if i == 0:
                ax1.legend(
                    facecolor="white",
                    fontsize=10,
                    framealpha=0,
                    bbox_to_anchor=(0, 1.3, 0.8, 0.2),
                    loc="upper center",
                    mode="expand",
                    borderaxespad=0,
                    ncol=1,
                )

        if save:
            self.save_plot(filename, bbox_inches="tight")

        return fig

    def plot_fom_dists(
        self,
        sigma_kerns: List[int],
        all_fom_dict: Dict[int, pd.Series],
        fom_limits: Dict[int, float],
        save: bool = False,
        filename: str = "all_fom",
    ):
        n = len(sigma_kerns)
        fig, axes = plt.subplots(n, constrained_layout=True)
        fig.set_figheight(n)
        fig.set_figwidth(5.5)
        fig.supylabel("Freq")

        xlim_lower = np.inf
        xlim_upper = -np.inf
        for sigma_kern, all_fom in all_fom_dict.items():
            xlim_lower = min(xlim_lower, min(all_fom))
            xlim_upper = max(xlim_upper, max(all_fom), fom_limits[sigma_kern])
        lims = PlotLimits(xlower=xlim_lower, xupper=xlim_upper)

        for i in range(n):
            if len(sigma_kerns) == 1:
                sigma_kern = sigma_kerns[0]
                ax: Axes = axes
            else:
                sigma_kern = sigma_kerns[i]
                ax: Axes = axes[i]
            all_fom = all_fom_dict[sigma_kern]
            fom_limit = fom_limits[sigma_kern]

            self._setup_ax(
                ax,
                lims,
                xlabel=False,
                xticks=i >= n - 1,
                ylabel=False,
                yticks=False,
            )
            if i >= n - 1:
                ax.set_xlabel(r"$\Sigma_{\rm FOM}$")
            # ax.set_ylabel("Freq")

            ax.text(
                0.02,
                0.95,
                r"$\sigma_{\rm kernel}$ = " + str(sigma_kern),
                ha="left",
                va="top",
                transform=ax.transAxes,
                fontsize=11,
            )
            ax.hist(
                all_fom,
                bins=np.linspace(min(all_fom), max(all_fom), 20),
                color=CONTROL_FOM_COLOR,
            )

            ax.axvline(fom_limit, linewidth=1.5, color="k", linestyle=FOM_LIMIT_LS)
            ax.text(
                fom_limit + 0.5,
                0.5 * ax.get_ylim()[1],
                r"$\Sigma_{\rm FOM, limit}$ = " + f"{fom_limit:0.2f}",
                fontsize=11,
            )

        if save:
            self.save_plot(filename, bbox_inches="tight")

        return fig


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
        control_index: bool = 0,
        custom_lims: Optional[PlotLimits] = None,
        title: str | None = None,
        save_filename: str = None,
    ):
        print(f"Plotting cut for flag {hex(flag)}...")
        fig = super().plot_cut(
            sn, flag, control_index, custom_lims, title, save_filename
        )
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
