#!/usr/bin/env python

from typing import Callable, List
import sys, argparse
import pandas as pd
import numpy as np
from copy import deepcopy
from lightcurve import (
    DEFAULT_CUT_NAMES,
    Cut,
    CutList,
    LimCutsTable,
    SnInfoTable,
    Supernova,
    AveragedSupernova,
    get_mjd0,
)
from download import (
    Credentials,
    load_config,
    make_dir_if_not_exists,
    parse_comma_separated_string,
)
from plot import PlotPdf

"""
UTILITY
"""


def hexstring_to_int(hexstring):
    return int(hexstring, 16)


class OutputReadMe:
    def __init__(self, output_dir, tnsname, cut_list, num_controls=0):
        filename = f"{output_dir}/{tnsname}/README.md"
        print(
            f"\nOpening README.md file for outputting cut information at {filename}..."
        )
        self.f = open(filename, "w+")
        self.tnsname: str = tnsname
        self.cut_list: CutList = cut_list
        self.begin(num_controls=num_controls)
        print("Success")

    def begin(self, num_controls=0):
        badday_cut = self.cut_list.get("badday_cut")
        if badday_cut:
            mjdbinsize = badday_cut.params["mjd_bin_size"]
        else:
            mjdbinsize = 1.0

        self.f.write(f"# SN {self.tnsname} Light Curve Cleaning and Averaging")
        self.f.write(
            f'\n\nThe ATLAS SN light curves are separated by filter (orange and cyan) and labelled as such in the file name. Averaged light curves contain an additional number in the file name that represents the MJD bin size used. Control light curves are located in the "controls" subdirectory and follow the same naming scheme, only with their control index added after the SN name.'
        )

        self.f.write(
            f"\n\nThe following details the file names for each of the light curve versions:"
        )
        self.f.write(
            f"\n\t- Original SN light curves: {self.tnsname}.o.lc.txt and {self.tnsname}.c.lc.txt"
        )
        self.f.write(
            f"\n\t- Cleaned SN light curves: {self.tnsname}.o.clean.lc.txt and {self.tnsname}.c.clean.lc.txt"
        )
        if self.cut_list.has("badday_cut"):
            self.f.write(
                f"\n\t- Averaged light curves (for MJD bin size {mjdbinsize:0.2f} days): {self.tnsname}.o.{mjdbinsize:0.2f}days.lc.txt and {self.tnsname}.c.{mjdbinsize:0.2f}days.lc.txt"
            )
        if self.cut_list.has("controls_cut"):
            self.f.write(
                f"\n\t- Control light curves, where X=001,...,{num_controls:03d}: {self.tnsname}_iX.o.lc.txt and {self.tnsname}_iX.c.lc.txt"
            )

        self.f.write(
            f'\n\nThe following summarizes the hex values in the "Mask" column of each light curve for each cut applied (see below sections for more information on each cut): '
        )
        if self.cut_list.has("uncert_cut"):
            self.f.write(
                f'\n\t- Uncertainty cut: {hex(self.cut_list.get("uncert_cut").flag)}'
            )
        if self.cut_list.has("x2_cut"):
            self.f.write(
                f'\n\t- Chi-square cut: {hex(self.cut_list.get("x2_cut").flag)}'
            )
        if self.cut_list.has("controls_cut"):
            self.f.write(
                f'\n\t- Control light curve cut: {hex(self.cut_list.get("controls_cut").flag)}'
            )
        if self.cut_list.has("badday_cut"):
            self.f.write(
                f'\n\t- Bad day (for averaged light curves only): {hex(self.cut_list.get("badday_cut").flag)}'
            )

        custom_cuts = self.cut_list.get_custom_cuts()
        for name in custom_cuts:
            self.f.write(
                f'\n\t- Custom cut {name[-1]} on "{self.cut_list.get(name).column}" column: {hex(self.cut_list.get(name).flag)}'
            )

    def add_filter_section(self, filt):
        self.f.write(f"\n\n## FILTER: {filt}")

    def add_template_correction_section(self):
        self.f.write(f"\n\n### ATLAS template change correction\n")
        # TODO

    def add_uncert_est_section(
        self,
        sigma_typical_old,
        sigma_typical_new,
        final_sigma_extra,
        percent_greater,
        apply,
    ):
        self.f.write(f"\n\n### True uncertainties estimation\n")
        self.f.write(
            f"\nWe can increase the typical uncertainties from {sigma_typical_old:0.2f} to {sigma_typical_new:0.2f} by adding an additional systematic uncertainty of {final_sigma_extra:0.2f} in quadrature."
        )
        self.f.write(
            f"\nNew typical uncertainty is {percent_greater:0.2f}% greater than old typical uncertainty."
        )
        self.f.write(f"\nApply true uncertainties estimation: {apply}.")
        if percent_greater >= 10:
            self.f.write("\nTrue uncertainties estimation recommended.")
            if apply:
                self.f.write(
                    '\nThe extra noise was added to the uncertainties of the SN light curve and copied to the "duJy_new" column.'
                )
        else:
            self.f.write(
                "\nTrue uncertainties estimation not needed; procedure skipped."
            )

    def add_uncert_cut_section(self, cut: Cut, percent_cut):
        self.f.write(f"\n\n### Uncertainty cut\n")
        self.f.write(
            f"\nTotal percent of SN light curve flagged with {hex(cut.flag)}: {percent_cut:0.2f}%"
        )

    def add_x2_cut_section(
        self, cut: Cut, percent_contamination, percent_loss, percent_cut
    ):
        self.f.write(f"\n\n### Chi-square cut\n")
        self.f.write(
            f"\nChi-square cut {cut.max_value:0.2f} selected with {percent_contamination:0.2f}% contamination and {percent_loss:0.2f}% loss\n"
        )
        self.f.write(
            f"\nTotal percent of SN light curve flagged with {hex(cut.flag)}: {percent_cut:0.2f}%"
        )

    def add_controls_cut_section(
        self,
        cut: Cut,
        x2_percent_cut,
        stn_percent_cut,
        Nclip_percent_cut,
        Ngood_percent_cut,
        questionable_percent_cut,
        percent_cut,
    ):
        self.f.write(f"\n\n### Control light curve cut\n")
        self.f.write(
            f'\nPercent of SN light curve above x2_max bound ({hex(cut.params["x2_flag"])}): {x2_percent_cut:0.2f}%'
        )
        self.f.write(
            f'\nPercent of SN light curve above stn_max bound ({hex(cut.params["stn_flag"])}): {stn_percent_cut:0.2f}%'
        )
        self.f.write(
            f'\nPercent of SN light curve above Nclip_max bound ({hex(cut.params["Nclip_flag"])}): {Nclip_percent_cut:0.2f}%'
        )
        self.f.write(
            f'\nPercent of SN light curve below Ngood_min bound ({hex(cut.params["Ngood_flag"])}): {Ngood_percent_cut:0.2f}%'
        )
        self.f.write(
            f'\nTotal percent of SN light curve flagged as questionable (not masked with control light curve flags but Nclip > 0) ({hex(cut.params["questionable_flag"])}): {questionable_percent_cut:0.2f}%'
        )
        self.f.write(
            f"\nTotal percent of SN light curve flagged as bad ({hex(cut.flag)}): {percent_cut:0.2f}%"
        )

    def add_badday_cut_section(self, percent_cut):
        self.f.write(
            f'\n\nAfter the cuts are applied, the light curves are resaved with the new "Mask" column.'
        )
        self.f.write(
            f"\nTotal percent of SN light curve flagged as bad ({hex(self.cut_list.get_all_flags())}): {percent_cut:0.2f}"
        )

        self.f.write(f"\n\n### Bad day cut (averaging)\n")
        # TODO

        self.f.write(
            f"\nThe averaged light curves are then saved in a new file with the MJD bin size added to the filename."
        )

    def add_custom_cut_section(self, name, cut: Cut, percent_cut):
        self.f.write(f"\n\n### Custom cut {name[-1]}\n")
        self.f.write(
            f"\nTotal percent of SN light curve flagged ({hex(cut.flag)}): {percent_cut:0.2f}%"
        )

    def save(self):
        self.f.close()


class UncertEstTable:
    def __init__(self, directory, filename=None):
        if filename is None:
            self.filename = f"{directory}/uncert_est_info.txt"
        else:
            self.filename = f"{directory}/{filename}"

        try:
            print(
                f"\nLoading true uncertainties estimation table at {self.filename}..."
            )
            self.t = pd.read_table(self.filename, delim_whitespace=True)
            print("Success")
        except:
            print(
                f"No existing true uncertainties estimation table; creating blank table..."
            )
            self.t = pd.DataFrame(
                columns=[
                    "tnsname",
                    "filter",
                    "sigma_extra",
                    "sigma_typical_old",
                    "sigma_typical_new",
                    "sigma_typical_new_pct_greater",
                    "recommended",
                    "applied",
                ]
            )

    def add_row(self, row):
        tnsname = row["tnsname"]
        filt = row["filter"]

        if len(self.t) > 0:
            matching_ix = np.where(
                self.t["tnsname"].eq(tnsname) & self.t["filter"].eq(filt)
            )[0]
            if len(matching_ix) > 1:
                raise RuntimeError(
                    f"ERROR: true uncertainties estimation table has {len(matching_ix)} matching rows for TNS name {tnsname} and filter {filt}"
                )

            if len(matching_ix) > 0:
                # update existing row
                idx = matching_ix[0]
                self.t.loc[idx, :] = row
            else:
                # new row
                self.t = pd.concat([self.t, pd.DataFrame([row])], ignore_index=True)
        else:
            # new row
            self.t = pd.concat([self.t, pd.DataFrame([row])], ignore_index=True)

    def save(self):
        print(f"\nSaving true uncertainties estimation table at {self.filename}...")
        self.t.to_string(self.filename)
        print("Success")


class ChiSquareCutTable:
    def __init__(self, directory, filename=None):
        if filename is None:
            self.filename = f"{directory}/x2_cut_info.txt"
        else:
            self.filename = f"{directory}/{filename}"

        try:
            print(f"\nLoading chi-square cut table at {self.filename}...")
            self.t = pd.read_table(self.filename, delim_whitespace=True)
            print("Success")
        except:
            print(f"No existing chi-square cut table; creating blank table...")
            self.t = pd.DataFrame(
                columns=[
                    "tnsname",
                    "filter",
                    "x2_cut",
                    "use_preSN_lc",
                    "stn_bound",
                    "pct_contamination",
                    "pct_loss",
                ]
            )

    def add_row(self, row):
        tnsname = row["tnsname"]
        filt = row["filter"]

        if len(self.t) > 0:
            matching_ix = np.where(
                self.t["tnsname"].eq(tnsname) & self.t["filter"].eq(filt)
            )[0]
            if len(matching_ix) > 1:
                raise RuntimeError(
                    f"ERROR: chi-square cut table has {len(matching_ix)} matching rows for TNS name {tnsname} and filter {filt}"
                )

            if len(matching_ix) > 0:
                # update existing row
                idx = matching_ix[0]
                self.t.loc[idx, :] = row
            else:
                # new row
                self.t = pd.concat([self.t, pd.DataFrame([row])], ignore_index=True)
        else:
            # new row
            self.t = pd.concat([self.t, pd.DataFrame([row])], ignore_index=True)

    def save(self):
        print(f"\nSaving chi-square cut table at {self.filename}...")
        self.t.to_string(self.filename)
        print("Success")


class CleanLoop:
    def __init__(
        self,
        input_dir,
        output_dir,
        credentials,
        sninfo_filename=None,
        flux2mag_sigmalimit=3.0,
        overwrite=False,
    ):
        self.sn: Supernova = None
        self.avg_sn: AveragedSupernova = None
        self.cut_list: CutList = None
        self.f: OutputReadMe = None
        self.p: PlotPdf = None

        self.credentials: Credentials = Credentials(
            credentials["atlas_username"],
            credentials["atlas_password"],
            credentials["tns_api_key"],
            credentials["tns_id"],
            credentials["tns_bot_name"],
        )
        self.input_dir: str = input_dir
        self.output_dir: str = output_dir
        self.flux2mag_sigmalimit: float = flux2mag_sigmalimit
        self.overwrite: bool = overwrite

        self.sninfo: SnInfoTable = SnInfoTable(
            self.output_dir, filename=sninfo_filename
        )
        self.uncert_est_info: UncertEstTable = UncertEstTable(self.output_dir)
        if cut_list.has("x2_cut"):
            self.x2_cut_info: ChiSquareCutTable = ChiSquareCutTable(self.output_dir)

    def apply_template_correction(self):
        print(f"\nApplying ATLAS template change correction:")
        # TODO: add_template_correction_section

    def check_uncert_est(self, cut: Cut, apply_function: Callable, plot: bool = False):
        print(f"\nChecking true uncertainties estimation:")

        stats = self.sn.get_uncert_est_stats(cut)
        final_sigma_extra = np.median(stats["sigma_extra"])

        sigma_typical_old = np.median(stats["median_dflux"])
        sigma_typical_new = np.sqrt(final_sigma_extra**2 + sigma_typical_old**2)
        percent_greater = 100 * (
            (sigma_typical_new - sigma_typical_old) / sigma_typical_old
        )
        print(
            f"We can increase the typical uncertainties from {sigma_typical_old:0.2f} to {sigma_typical_new:0.2f} by adding an additional systematic uncertainty of {final_sigma_extra:0.2f} in quadrature"
        )
        print(
            f"New typical uncertainty is {percent_greater:0.2f}% greater than old typical uncertainty"
        )

        apply = apply_function()
        print(f"Apply true uncertainties estimation: {apply}")
        if percent_greater >= 10:
            print("True uncertainties estimation recommended")
            print(f'{"Applying" if apply else "Skipping"} procedure...')
            if apply:
                self.sn.add_noise_to_dflux(final_sigma_extra)
                print("Success")
                print(
                    'The extra noise was added to the uncertainties of the SN light curve and copied to the "duJy_new" column'
                )

                if plot:
                    self.p.plot_uncert_est(
                        self.sn.lcs[0],
                        self.sn.tnsname,
                        self.p.get_lims(
                            lc=self.sn.lcs[0],
                            indices=self.sn.lcs[0].get_good_indices(
                                self.cut_list.get("uncert_cut").flag
                            ),
                        ),
                    )
        else:
            print("True uncertainties estimation not needed; skipping procedure...")

        self.f.add_uncert_est_section(
            sigma_typical_old,
            sigma_typical_new,
            final_sigma_extra,
            percent_greater,
            apply,
        )

        uncert_est_info_row = {
            "tnsname": self.sn.tnsname,
            "filter": self.sn.filt,
            "sigma_extra": final_sigma_extra,
            "sigma_typical_old": sigma_typical_old,
            "sigma_typical_new": sigma_typical_new,
            "sigma_typical_new_pct_greater": percent_greater,
            "recommended": percent_greater >= 10,
            "applied": apply,
        }
        return apply, uncert_est_info_row

    def apply_uncert_cut(self, cut: Cut, plot: bool = False):
        if cut is None:
            return
        print(f"\nApplying uncertainty cut ({cut}):")
        percent_cut = self.sn.apply_cut(cut)
        print("Success")
        print(
            f"Total percent of SN light curve flagged with {hex(cut.flag)}: {percent_cut:0.2f}%"
        )
        self.f.add_uncert_cut_section(cut, percent_cut)

        if plot:
            self.p.plot_cut(
                self.sn.lcs[0],
                cut.flag,
                self.p.get_lims(
                    lc=self.sn.lcs[0],
                    indices=self.sn.lcs[0].get_good_indices(cut.flag),
                ),
                title="Uncertainty cut",
            )

    def apply_x2_cut(self, cut: Cut, plot: bool = False):
        if cut is None:
            return
        print(f"\nApplying chi-square cut ({cut}):")

        if cut.params["use_pre_mjd0_lc"]:
            print("Using pre-MJD0 light curve to determine contamination and loss")
            if self.sn.mjd0 is None:
                raise RuntimeError(
                    "ERROR: MJD0 cannot be None. Please provide MJD0 throught the SN info table or --mjd0 argument, or set the use_pre_MJD0_lc field in the config file to False."
                )
            lc_temp = deepcopy(self.sn.lcs[0])
            ix = lc_temp.ix_inrange("MJD", uplim=self.sn.mjd0)
        else:
            print("Using control light curves to determine contamination and loss")
            if self.sn.num_controls < 1:
                raise RuntimeError(
                    "ERROR: No control light curves loaded. Use the --num_controls argument to load control light curves, or change the [x2_cut][use_pre_mjd0_lc] field to True."
                )
            lc_temp = self.sn.get_all_controls()
            ix = lc_temp.t.index.values

        limcuts = LimCutsTable(lc_temp, cut.params["stn_bound"], indices=ix)
        limcuts.calculate_table(
            cut.params["min_cut"], cut.params["max_cut"], cut.params["cut_step"]
        )
        print("Success")

        data = limcuts.calculate_row(cut.max_value)
        print(
            f'Applying chi-square cut of {cut.max_value:0.2f} with {data["Pcontamination"]:0.2f}% contamination and {data["Ploss"]:0.2f}% loss...'
        )
        percent_cut = self.sn.apply_cut(cut)
        print("Success")
        print(
            f"Total percent of SN light curve flagged with {hex(cut.flag)}: {percent_cut:0.2f}%"
        )

        if plot:
            self.p.plot_cut(
                self.sn.lcs[0],
                cut.flag,
                self.p.get_lims(
                    lc=self.sn.lcs[0],
                    indices=self.sn.lcs[0].get_good_indices(cut.flag),
                ),
                title="Chi-square cut",
            )

        self.f.add_x2_cut_section(
            cut, data["Pcontamination"], data["Ploss"], percent_cut
        )

        x2_info_row = {
            "tnsname": self.sn.tnsname,
            "filter": self.sn.filt,
            "x2_cut": cut.max_value,
            "use_pre_mjd0_lc": cut.params["use_pre_mjd0_lc"],
            "stn_bound": cut.params["stn_bound"],
            "pct_contamination": round(data["Pcontamination"], 2),
            "pct_loss": round(data["Ploss"], 2),
        }
        return x2_info_row

    def apply_controls_cut(self, cut: Cut, previous_flags: int, plot: bool = False):
        if cut is None:
            return
        print(f"\nApplying control light curve cut ({cut}):")

        (
            x2_percent_cut,
            stn_percent_cut,
            Nclip_percent_cut,
            Ngood_percent_cut,
            questionable_percent_cut,
            percent_cut,
        ) = self.sn.apply_controls_cut(cut, previous_flags)
        print("Success")

        print(
            f'Percent of data above x2_max bound ({hex(cut.params["x2_flag"])}): {x2_percent_cut:0.2f}%'
        )
        print(
            f'Percent of data above stn_max bound ({hex(cut.params["stn_flag"])}): {stn_percent_cut:0.2f}%'
        )
        print(
            f'Percent of data above Nclip_max bound ({hex(cut.params["Nclip_flag"])}): {Nclip_percent_cut:0.2f}%'
        )
        print(
            f'Percent of data below Ngood_min bound ({hex(cut.params["Ngood_flag"])}): {Ngood_percent_cut:0.2f}%'
        )
        print(
            f'Total percent of data flagged as questionable (not masked with control light curve flags but Nclip > 0) ({hex(cut.params["questionable_flag"])}): {questionable_percent_cut:0.2f}%'
        )
        print(
            f"Total percent of data flagged as bad ({hex(cut.flag)}): {percent_cut:0.2f}%"
        )
        self.f.add_controls_cut_section(
            cut,
            x2_percent_cut,
            stn_percent_cut,
            Nclip_percent_cut,
            Ngood_percent_cut,
            questionable_percent_cut,
            percent_cut,
        )

        if plot:
            self.p.plot_cut(
                self.sn.lcs[0],
                cut.flag,
                self.p.get_lims(
                    lc=self.sn.lcs[0],
                    indices=self.sn.lcs[0].get_good_indices(cut.flag),
                ),
                title="Control light curve cut",
            )

    def apply_badday_cut(self, cut: Cut, previous_flags, plot: bool = False):
        if cut is None:
            return
        print(f"\nApplying bad day cut (averaging) ({cut}):")
        self.avg_sn, percent_cut = self.sn.apply_badday_cut(
            cut, previous_flags, flux2mag_sigmalimit=self.flux2mag_sigmalimit
        )
        print("Success")

        print(
            f"Total percent of SN light curve flagged as bad ({hex(self.cut_list.get_all_flags())}): {percent_cut:0.2f}"
        )
        self.f.add_badday_cut_section(percent_cut)

        if plot:
            lims = self.p.get_lims(
                lc=self.avg_sn.avg_lcs[0],
                indices=self.avg_sn.avg_lcs[0].get_good_indices(cut.flag),
            )

            self.p.plot_cut(
                self.avg_sn.avg_lcs[0],
                cut.flag,
                lims,
                title="Bad day cut",
            )

            self.p.plot_averaged_SN(
                self.avg_sn, cut.flag, lims, plot_controls=True, plot_flagged=False
            )

    def apply_custom_cut(self, name, cut: Cut, plot: bool = False):
        cut = self.cut_list.get(name)
        print(f"\nApplying custom cut ({cut})...")
        percent_cut = self.sn.apply_cut(cut)
        print("Success")

        print(
            f"Total percent of SN light curve flagged with {hex(cut.flag)}: {percent_cut:0.2f}%"
        )
        self.f.add_custom_cut_section(name, cut, percent_cut)

        if plot:
            self.p.plot_cut(
                self.sn.lcs[0],
                cut.flag,
                self.p.get_lims(
                    lc=self.sn.lcs[0],
                    indices=self.sn.lcs[0].get_good_indices(cut.flag),
                ),
                title=f"Custom cut {name}",
            )

    def clean_lcs(
        self,
        tnsname: str,
        mjd0,
        filt: str,
        apply_uncert_est_function: Callable,
        num_controls: int = 0,
        apply_template_correction: bool = False,
        plot: bool = False,
    ):
        print(f"\n\tFILTER: {filt}")

        # load the SN and control light curves
        self.sn = Supernova(tnsname=tnsname, mjd0=mjd0, filt=filt)
        try:
            self.sn.load_all(self.input_dir, num_controls=num_controls)
        except Exception as e:
            raise RuntimeError(f"ERROR: Could not load light curves: {str(e)}")

        # prepare the light curves for cleaning
        print()
        self.sn.prep_for_cleaning(verbose=True)

        if plot:
            # initialize PDF of diagnostic plots
            self.p = PlotPdf(f"{self.output_dir}/{tnsname}", tnsname, filt=filt)
            self.p.plot_SN(
                self.sn,
                self.p.get_lims(lc=self.sn.lcs[0]),
                plot_controls=True,
                plot_template_changes=True,
            )

        # template correction
        if apply_template_correction:
            self.apply_template_correction()

        # uncertainty cut
        self.apply_uncert_cut(self.cut_list.get("uncert_cut"), plot=plot)

        # true uncertainties estimation
        _, uncert_est_info_row = self.check_uncert_est(
            self.cut_list.get("uncert_est"),
            apply_function=apply_uncert_est_function,
            plot=plot,
        )
        self.uncert_est_info.add_row(uncert_est_info_row)

        # chi-square cut
        x2_info_row = self.apply_x2_cut(self.cut_list.get("x2_cut"), plot=plot)
        if self.cut_list.has("x2_cut"):
            self.x2_cut_info.add_row(x2_info_row)

        # control light curve cut
        self.apply_controls_cut(
            self.cut_list.get("controls_cut"),
            previous_flags=self.cut_list.get_previous_flags("controls_cut"),
            plot=plot,
        )

        # custom cuts
        custom_cuts = self.cut_list.get_custom_cuts()
        for name, cut in custom_cuts.items():
            self.apply_custom_cut(name, cut, plot=plot)

        # plot the cleaned light curves so far
        previous_flags = self.cut_list.get_previous_flags("badday_cut")
        if plot:
            lims = self.p.get_lims(
                lc=self.sn.lcs[0],
                indices=self.sn.lcs[0].get_good_indices(previous_flags),
            )
            self.p.plot_cut(
                self.sn.lcs[0], previous_flags, lims, title="All previous cuts"
            )
            self.p.plot_cleaned_SN(
                self.sn, previous_flags, lims, plot_controls=True, plot_flagged=False
            )

        # bad day cut (averaging)
        self.apply_badday_cut(
            self.cut_list.get("badday_cut"),
            previous_flags=previous_flags,
            plot=plot,
        )

        # save cleaned SN and control light curves
        self.sn.save_all(self.output_dir, overwrite=self.overwrite)

        if self.cut_list.has("badday_cut"):
            # save averaged SN and control light curves
            self.avg_sn.save_all(self.output_dir, overwrite=self.overwrite)

        if self.cut_list.has("uncert_est"):
            # save uncertainty estimation table
            self.uncert_est_info.save()

        if self.cut_list.has("x2_cut"):
            # save chi-square cut table
            self.x2_cut_info.save()

        # save the SN info table
        self.sninfo.save()

        if plot:
            # save the PDF of diagnostic plots
            self.p.save_pdf()

    def loop(
        self,
        tnsnames: List[str],
        apply_uncert_est_function: Callable,
        num_controls: int = 0,
        mjd0=None,
        filters: List[str] = ["o", "c"],
        cut_list: CutList = None,
        apply_template_correction: bool = False,
        plot: bool = False,
    ):
        self.cut_list = cut_list

        for obj_index in range(len(tnsnames)):
            tnsname = tnsnames[obj_index]
            print(f"\n\tCLEANING LIGHT CURVES FOR: SN {tnsname}")

            make_dir_if_not_exists(f"{output_dir}/{tnsname}")
            self.f = OutputReadMe(
                self.output_dir, tnsname, cut_list, num_controls=num_controls
            )

            if mjd0 is None and (
                plot or cut_list.get("x2_cut").params["use_pre_mjd0_lc"]
            ):
                mjd0, coords = get_mjd0(tnsname, self.sninfo, self.credentials)
                if not coords is None:
                    print(f"Setting MJD0 to {mjd0}")
                    self.sninfo.update_row(tnsname, coords=coords, mjd0=mjd0)
                    print("Success")
            else:
                print(f"\nSetting MJD0 to {mjd0}")

            for filt in filters:
                self.f.add_filter_section(filt)
                self.clean_lcs(
                    tnsname,
                    mjd0,
                    filt,
                    apply_uncert_est_function,
                    num_controls=num_controls,
                    apply_template_correction=apply_template_correction,
                    plot=plot,
                )


def parse_config_filters(args, config):
    if args.filters:
        return parse_comma_separated_string(args.filters)
    else:
        return parse_comma_separated_string(config["convert"]["filters"])


def find_config_custom_cuts(config):
    print("\nSearching config file for custom cuts...")
    custom_cuts = []
    for key in config:
        if key.endswith("_cut") and not key in DEFAULT_CUT_NAMES:
            custom_cuts.append(config[key])
    print(f"Found {len(custom_cuts)}")
    return custom_cuts


def parse_config_cuts(args, config):
    cut_list = CutList()
    if args.custom_cuts:
        config_custom_cuts = find_config_custom_cuts(config)

    print(f"\nProcedures to apply:")

    # always check true uncertainties estimation, but will only apply if args.true_uncert_est
    temp_x2_max_value = float(config["uncert_est"]["temp_x2_max_value"])
    print(
        f"- True uncertainties estimation check: temporary chi-square cut at {temp_x2_max_value}"
    )
    params = {
        "temp_x2_max_value": temp_x2_max_value,
        "uncert_cut_flag": hexstring_to_int(config["uncert_cut"]["flag"]),
    }
    uncert_est = Cut(params=params)
    cut_list.add(uncert_est, "uncert_est")

    if args.uncert_cut:
        uncert_cut = Cut(
            column="duJy",
            max_value=float(config["uncert_cut"]["max_value"]),
            flag=hexstring_to_int(config["uncert_cut"]["flag"]),
        )
        cut_list.add(uncert_cut, "uncert_cut")
        print(f"- Uncertainty cut: {uncert_cut}")

    if args.x2_cut:
        params = {
            "stn_bound": float(config["x2_cut"]["stn_bound"]),
            "min_cut": int(config["x2_cut"]["min_cut"]),
            "max_cut": int(config["x2_cut"]["max_cut"]),
            "cut_step": int(config["x2_cut"]["cut_step"]),
            "use_pre_mjd0_lc": config["x2_cut"]["use_pre_mjd0_lc"] == "True",
        }
        x2_cut = Cut(
            column="chi/N",
            max_value=float(config["x2_cut"]["max_value"]),
            flag=hexstring_to_int(config["x2_cut"]["flag"]),
            params=params,
        )
        cut_list.add(x2_cut, "x2_cut")
        print(f"- Chi-square cut: {x2_cut}")

    if args.controls_cut:
        params = {
            "questionable_flag": hexstring_to_int(
                config["controls_cut"]["questionable_flag"]
            ),
            "x2_max": float(config["controls_cut"]["x2_max"]),
            "x2_flag": hexstring_to_int(config["controls_cut"]["x2_flag"]),
            "stn_max": float(config["controls_cut"]["stn_max"]),
            "stn_flag": hexstring_to_int(config["controls_cut"]["stn_flag"]),
            "Nclip_max": int(config["controls_cut"]["Nclip_max"]),
            "Nclip_flag": hexstring_to_int(config["controls_cut"]["Nclip_flag"]),
            "Ngood_min": int(config["controls_cut"]["Ngood_min"]),
            "Ngood_flag": hexstring_to_int(config["controls_cut"]["Ngood_flag"]),
        }
        controls_cut = Cut(
            flag=hexstring_to_int(config["controls_cut"]["bad_flag"]), params=params
        )
        cut_list.add(controls_cut, "controls_cut")
        print(f"- Control light curve cut: {controls_cut}")

    if args.averaging:
        params = {
            "mjd_bin_size": (
                float(config["averaging"]["mjd_bin_size"])
                if args.mjd_bin_size is None
                else args.args.mjd_bin_size
            ),
            "x2_max": float(config["averaging"]["x2_max"]),
            "Nclip_max": int(config["averaging"]["Nclip_max"]),
            "Ngood_min": int(config["averaging"]["Ngood_min"]),
            "ixclip_flag": hexstring_to_int(config["averaging"]["ixclip_flag"]),
            "smallnum_flag": hexstring_to_int(config["averaging"]["smallnum_flag"]),
        }
        badday_cut = Cut(
            flag=hexstring_to_int(config["averaging"]["flag"]), params=params
        )
        cut_list.add(badday_cut, "badday_cut")
        print(f"- Bad day cut (averaging): {badday_cut}")

    if args.custom_cuts:
        for i in range(len(config_custom_cuts)):
            cut_settings = config_custom_cuts[i]
            try:
                custom_cut = Cut(
                    column=cut_settings["column"],
                    flag=hexstring_to_int(cut_settings["flag"]),
                    min_value=(
                        float(cut_settings["min_value"])
                        if cut_settings["min_value"] != "None"
                        else None
                    ),
                    max_value=(
                        float(cut_settings["max_value"])
                        if cut_settings["max_value"] != "None"
                        else None
                    ),
                )
                cut_list.add(custom_cut, f"custom_cut_{i}")
                print(f"- Custom cut {i}: {custom_cut}")
            except Exception as e:
                print(f"WARNING: Could not parse custom cut {cut_settings}: {str(e)}")

    has_duplicate_flags, duplicate_flags = cut_list.check_for_flag_duplicates()
    if has_duplicate_flags:
        raise RuntimeError(
            f"ERROR: Cuts in the config file contain duplicate flags: {duplicate_flags}."
        )
    return cut_list


# define command line arguments
def define_args(parser=None, usage=None, conflict_handler="resolve"):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

    parser.add_argument(
        "tnsnames", nargs="+", help="TNS names of the transients to clean"
    )
    parser.add_argument(
        "--sninfo_file",
        default=None,
        type=str,
        help="file name of .txt file with SN info table",
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
        help="comma-separated list of filters to clean",
    )
    parser.add_argument(
        "-p",
        "--plot",
        default=False,
        action="store_true",
        help="store a summary PDF file of diagnostic plots",
    )

    # cleaning a single SN and/or controls
    parser.add_argument(
        "--mjd0", type=float, default=None, help="transient start date in MJD"
    )

    # cleaning control light curves
    # parser.add_argument('-c','--controls', default=False, action='store_true', help='clean control light curves in addition to transient light curve')
    parser.add_argument(
        "--num_controls",
        type=int,
        default=None,
        help="number of control light curves to load and clean",
    )

    # possible cuts
    parser.add_argument(
        "-t",
        "--template_correction",
        default=False,
        action="store_true",
        help="apply automatic ATLAS template change correction",
    )
    parser.add_argument(
        "-e",
        "--uncert_est",
        default=False,
        action="store_true",
        help="apply true uncertainty estimation",
    )
    parser.add_argument(
        "-u",
        "--uncert_cut",
        default=False,
        action="store_true",
        help="apply uncertainty cut",
    )
    parser.add_argument(
        "-x",
        "--x2_cut",
        default=False,
        action="store_true",
        help="apply chi-square cut",
    )
    parser.add_argument(
        "-c",
        "--controls_cut",
        default=False,
        action="store_true",
        help="apply control light curve cut",
    )
    parser.add_argument(
        "-g",
        "--averaging",
        default=False,
        action="store_true",
        help="average light curves and cut bad days",
    )
    parser.add_argument(
        "-m",
        "--mjd_bin_size",
        type=float,
        default=None,
        help="MJD bin size in days for averaging",
    )
    parser.add_argument(
        "--custom_cuts",
        default=False,
        action="store_true",
        help="scan config file for custom cuts",
    )

    return parser


if __name__ == "__main__":
    args = define_args().parse_args()
    config = load_config(args.config_file)

    if len(args.tnsnames) < 1:
        raise RuntimeError("ERROR: Please specify at least one TNS name to clean.")
    if len(args.tnsnames) > 1 and not args.mjd0 is None:
        raise RuntimeError(
            f"ERROR: Cannot specify one MJD0 {args.mjd0} for a batch of SNe."
        )
    print(f"\nList of transients to clean: {args.tnsnames}")

    input_dir = config["dir"]["atclean_input"]
    output_dir = config["dir"]["output"]
    sninfo_filename = config["dir"]["sninfo_filename"]
    make_dir_if_not_exists(input_dir)
    make_dir_if_not_exists(output_dir)
    print(f"\nATClean input directory: {input_dir}")
    print(f"Output directory: {output_dir}")

    print(f'TNS ID: {config["credentials"]["tns_id"]}')
    print(f'TNS bot name: {config["credentials"]["tns_bot_name"]}')

    print(f"Overwrite existing files: {args.overwrite}")
    print(f"Save PDF of diagnostic plots: {args.plot}")
    filters = parse_config_filters(args, config)
    print(f"Filters: {filters}")
    flux2mag_sigmalimit = float(config["download"]["flux2mag_sigmalimit"])
    print(f"Sigma limit when converting flux to magnitude: {flux2mag_sigmalimit}")
    if args.mjd0:
        print(f"MJD0: {args.mjd0}")

    # print(f'\nApplyin control light curve cut: {args.controls}')
    num_controls = (
        args.num_controls
        if args.num_controls
        else int(config["download"]["num_controls"])
    )
    print(f"Number of control light curves to clean: {num_controls}")

    cut_list = parse_config_cuts(args, config)

    print()
    clean = CleanLoop(
        input_dir,
        output_dir,
        config["credentials"],
        sninfo_filename=sninfo_filename,
        flux2mag_sigmalimit=flux2mag_sigmalimit,
        overwrite=args.overwrite,
    )

    def apply_uncert_est_function():
        return args.uncert_est

    clean.loop(
        args.tnsnames,
        apply_uncert_est_function,
        cut_list=cut_list,
        num_controls=num_controls,
        mjd0=args.mjd0,
        filters=filters,
        apply_template_correction=args.template_correction,
        plot=args.plot,
    )
