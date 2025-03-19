#!/usr/bin/env python

import os
from typing import List
import matplotlib
from matplotlib import gridspec
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from lightcurve import (
    TEMPLATE_CHANGE_1_MJD,
    TEMPLATE_CHANGE_2_MJD,
    LimCutsTable,
    Cut,
    LightCurve,
    Supernova,
    AveragedSupernova,
)

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

    def calc_ylims(
        self, lc: LightCurve | None = None, indices: List[int] | None = None
    ):
        if lc is None:
            print("No light curve provided; skipping plot limits calculation...")
            return

        if indices is None or len(indices) < 2:
            indices = lc.getindices()

        flux_min = lc.t.loc[indices, lc.colnames.flux].min()
        flux_max = lc.t.loc[indices, lc.colnames.flux].max()
        offset = 0.05 * abs(flux_max - flux_min)

        if self.ylower is None:
            self.ylower = flux_min - offset
        if self.yupper is None:
            self.yupper = flux_max + offset

    def get_xlims(self):
        return [self.xlower, self.xupper]

    def get_ylims(self):
        return [self.ylower, self.yupper]

    def is_empty(self):
        return (
            self.xlower is None
            and self.xupper is None
            and self.ylower is None
            and self.yupper is None
        )

    def __str__(self):
        return f"Plot limits: x-axis [{self.xlower}, {self.xupper}], y-axis [{self.ylower}, {self.yupper}]"


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
        lc: LightCurve = None,
        indices: List[int] = None,
        custom_lims: PlotLimits | None = None,
    ) -> PlotLimits:
        if custom_lims is not None:
            lims = PlotLimits(
                xlower=custom_lims.xlower,
                xupper=custom_lims.xupper,
                ylower=custom_lims.ylower,
                yupper=custom_lims.yupper,
            )
        else:
            lims = PlotLimits()

        lims.calc_ylims(lc=lc, indices=indices)
        return lims

    def plot_SN(
        self,
        sn: Supernova,
        lims: PlotLimits,
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
        ax1.set_xlabel(sn.colnames_master.mjd)
        ax1.axhline(linewidth=1, color="k")

        if plot_controls and sn.num_controls > 0:
            # plot control light curves
            label = f"{sn.num_controls} control light curves"
            for control_index in sn.get_control_indices():
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

        ax1.set_xlim(lims.xlower, lims.xupper)
        ax1.set_ylim(lims.ylower, lims.yupper)
        ax1.legend(loc="upper right", facecolor="white", framealpha=1.0).set_zorder(100)

        if save:
            self.save_plot(filename)

        return fig

    def plot_cut(
        self,
        lc: LightCurve,
        flag: int,
        lims: PlotLimits,
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
        ax1.set_xlabel(lc.colnames.mjd)
        ax2.axhline(linewidth=1, color="k")

        good_ix = lc.get_good_indices(flag)
        bad_ix = lc.get_bad_indices(flag)

        if lc.can_plot(good_ix):
            ax1.errorbar(
                lc.t.loc[good_ix, lc.colnames.mjd],
                lc.t.loc[good_ix, lc.colnames.flux],
                yerr=lc.t.loc[good_ix, lc.colnames.dflux_new],
                fmt="none",
                ecolor=SN_FLUX_COLORS[lc.filt],
                elinewidth=1,
                capsize=1.2,
                c=SN_FLUX_COLORS[lc.filt],
                alpha=0.5,
            )
            ax1.scatter(
                lc.t.loc[good_ix, lc.colnames.mjd],
                lc.t.loc[good_ix, lc.colnames.flux],
                s=marker_size,
                lw=marker_edgewidth,
                color=SN_FLUX_COLORS[lc.filt],
                marker="o",
                alpha=0.5,
                label="Cleaned measurements",
            )

            ax2.errorbar(
                lc.t.loc[good_ix, lc.colnames.mjd],
                lc.t.loc[good_ix, lc.colnames.flux],
                yerr=lc.t.loc[good_ix, lc.colnames.dflux_new],
                fmt="none",
                ecolor=SN_FLUX_COLORS[lc.filt],
                elinewidth=1,
                capsize=1.2,
                c=SN_FLUX_COLORS[lc.filt],
                alpha=0.5,
                zorder=5,
            )
            ax2.scatter(
                lc.t.loc[good_ix, lc.colnames.mjd],
                lc.t.loc[good_ix, lc.colnames.flux],
                s=marker_size,
                lw=marker_edgewidth,
                color=SN_FLUX_COLORS[lc.filt],
                marker="o",
                alpha=0.5,
                label="Cleaned measurements",
                zorder=5,
            )

        if lc.can_plot(bad_ix):
            ax2.errorbar(
                lc.t.loc[bad_ix, lc.colnames.mjd],
                lc.t.loc[bad_ix, lc.colnames.flux],
                yerr=lc.t.loc[bad_ix, lc.colnames.dflux_new],
                fmt="none",
                ecolor=SN_FLAGGED_FLUX_COLOR,
                elinewidth=1,
                capsize=1.2,
                c=SN_FLAGGED_FLUX_COLOR,
                alpha=0.5,
                zorder=10,
            )
            ax2.scatter(
                lc.t.loc[bad_ix, lc.colnames.mjd],
                lc.t.loc[bad_ix, lc.colnames.flux],
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

        ax1.set_xlim(lims.xlower, lims.xupper)
        ax1.set_ylim(lims.ylower, lims.yupper)
        ax2.set_xlim(lims.xlower, lims.xupper)
        ax2.set_ylim(lims.ylower, lims.yupper)

        ax1.legend(loc="upper right", facecolor="white", framealpha=1.0).set_zorder(100)
        ax2.legend(loc="upper right", facecolor="white", framealpha=1.0).set_zorder(100)

        if not save_filename is None:
            self.save_plot(save_filename)

        return fig

    def plot_cleaned_SN(
        self,
        sn: Supernova,
        flag: int,
        lims: PlotLimits,
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
        ax1.set_xlabel(sn.colnames_master.mjd)
        ax1.axhline(linewidth=1, color="k")

        if plot_controls and sn.num_controls > 0:
            # plot control light curves
            label = f"Cleaned control measurements"
            for control_index in sn.get_control_indices():
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

        ax1.set_xlim(lims.xlower, lims.xupper)
        ax1.set_ylim(lims.ylower, lims.yupper)
        ax1.legend(loc="upper right", facecolor="white", framealpha=1.0).set_zorder(100)

        if save:
            self.save_plot(filename)

        return fig

    def plot_averaged_SN(
        self,
        avg_sn: AveragedSupernova,
        flag: int,
        lims: PlotLimits,
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
        ax1.set_xlabel(avg_sn.colnames_master.mjd)
        ax1.axhline(linewidth=1, color="k")

        if plot_controls and avg_sn.num_controls > 0:
            # plot control light curves
            label = f"Cleaned & averaged control measurements"
            for control_index in avg_sn.get_control_indices():
                lc = avg_sn.avg_lcs[control_index]
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

        avg_sn_lc = avg_sn.avg_lcs[0]
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

        ax1.set_xlim(lims.xlower, lims.xupper)
        ax1.set_ylim(lims.ylower, lims.yupper)
        ax1.legend(loc="upper right", facecolor="white", framealpha=1.0).set_zorder(100)

        if save:
            self.save_plot(filename)

        return fig

    def plot_limcuts(
        self,
        limcuts: LimCutsTable,
        cut: Cut,
        cut_start: int,
        cut_stop: int,
        use_preSN_lc=False,
    ):
        # TODO
        pass

    def plot_uncert_est(
        self,
        lc: LightCurve,
        tnsname: str,
        lims: PlotLimits,
        save: bool = False,
        filename: str = "uncert_est",
    ):
        if not f"{lc.colnames.dflux}_new" in lc.t.columns:
            print(
                f"WARNING: Cannot plot true uncertainties estimation due to missing {lc.colnames.dflux}_new column; skipping..."
            )
            return None

        fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True)
        fig.set_figwidth(7)
        fig.set_figheight(5)

        ax1.set_title(
            f"SN {tnsname} {lc.filt}-band flux\nbefore true uncertainties estimation"
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

        ax1.set_xlim(lims.xlower, lims.xupper)
        ax1.set_ylim(lims.ylower, lims.yupper)
        ax2.set_xlim(lims.xlower, lims.xupper)
        ax2.set_ylim(lims.ylower, lims.yupper)

        if save:
            self.save_plot(filename)

        return fig

    def plot_template_correction(
        self,
        lc: LightCurve,
        lims: PlotLimits,
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
        ax1.set_xlim(lims.xlower, lims.xupper)
        ax1.set_ylim(lims.ylower, lims.yupper)

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
        ax2.set_ylim(lims.ylower, lims.yupper)

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
        ax3.set_ylim(lims.ylower, lims.yupper)

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
        lims: PlotLimits,
        plot_controls: bool = True,
        plot_template_changes: bool = True,
        save: bool = False,
        filename: str = "original",
    ):
        print(
            f'Plotting original SN{" and control light curves" if plot_controls else ""}...'
        )
        fig = super().plot_SN(
            sn, lims, plot_controls, plot_template_changes, save, filename
        )
        self.pdf.savefig(fig)

    def plot_cut(
        self,
        lc: LightCurve,
        flag: int,
        lims: PlotLimits,
        title: str | None = None,
        save_filename: str = None,
    ):
        print(f"Plotting cut for flag {hex(flag)}...")
        fig = super().plot_cut(lc, flag, lims, title, save_filename)
        self.pdf.savefig(fig)

    def plot_cleaned_SN(
        self,
        sn: Supernova,
        flag: int,
        lims: PlotLimits,
        plot_controls: bool = True,
        plot_flagged: bool = True,
        save: bool = False,
        filename: str = "cleaned",
    ):
        print(
            f'Plotting cleaned SN{" and control light curves" if plot_controls else ""} using flag {hex(flag)}...'
        )
        fig = super().plot_cleaned_SN(
            sn, flag, lims, plot_controls, plot_flagged, save, filename
        )
        self.pdf.savefig(fig)

    def plot_averaged_SN(
        self,
        avg_sn: AveragedSupernova,
        flag: int,
        lims: PlotLimits,
        plot_controls: bool = True,
        plot_flagged: bool = True,
        save: bool = False,
        filename: str = "averaged",
    ):
        print(
            f'Plotting averaged SN{" and control light curves" if plot_controls else ""} using flag {hex(flag)}...'
        )
        fig = super().plot_averaged_SN(
            avg_sn, flag, lims, plot_controls, plot_flagged, save, filename
        )
        self.pdf.savefig(fig)

    def plot_limcuts(
        self,
        limcuts: LimCutsTable,
        cut: Cut,
        cut_start: int,
        cut_stop: int,
        use_preSN_lc=False,
    ):
        print("Plotting LimCutsTable...")
        fig = super().plot_limcuts(limcuts, cut, cut_start, cut_stop, use_preSN_lc)
        self.pdf.savefig(fig)

    def plot_uncert_est(
        self,
        lc: LightCurve,
        tnsname: str,
        lims: PlotLimits,
        save: bool = False,
        filename: str = "uncert_est",
    ):
        print("Plotting true uncertainties estimation...")
        fig = super().plot_uncert_est(lc, tnsname, lims, save, filename)
        if not fig is None:
            self.pdf.savefig(fig)

    def plot_template_correction(self, lc: LightCurve):
        print("Plotting ATLAS template chanages correction...")
        fig = super().plot_template_correction(lc)
        self.pdf.savefig(fig)
