#!/usr/bin/env python

from datetime import datetime
from typing import Callable, Dict, List, Optional
import sys, argparse
import pandas as pd
import numpy as np
from copy import deepcopy
from lightcurve import (
    LimCutsTable,
    Supernova,
    AveragedSupernova,
)
from download import (
    Credentials,
)
from plot import PlotPdf
from utils import (
    CONFIG_CUT_NAMES,
    BadDayCut,
    ChiSquareCut,
    ControlLightCurveCut,
    CustomCut,
    UncertaintyCut,
    UncertaintyEstimation,
    find_all_filts,
    format_float_string,
    get_config_custom_cuts,
    hexstring_to_int,
    Cut,
    CutList,
    SnInfoTable,
    get_allowed_presets,
    get_mjd0_from_tns,
    PresetColumnNames,
    new_row,
    parse_config_str,
    load_config,
    make_dir_if_not_exists,
    parse_comma_separated_string,
)


class OutputReadMe:
    def __init__(self, output_dir, tnsname, cut_list, num_controls=0):
        timestamp = datetime.now()
        filename = (
            f"{output_dir}/{tnsname}/README_{timestamp.strftime('%Y%m%d_%H%M%S')}.md"
        )
        print(
            f"\nOpening README.md file for outputting cut information at {filename}..."
        )
        self.f = open(filename, "w+")
        self.tnsname: str = tnsname
        self.cut_list: CutList = cut_list
        self.begin(timestamp, num_controls=num_controls)
        print("Success")

    def begin(self, timestamp: datetime, num_controls: int = 0):
        badday_cut = self.cut_list.get(BadDayCut.name())
        mjdbinsize = 1.0 if badday_cut is None else badday_cut.mjd_bin_size

        self.f.write(f"# SN {self.tnsname} Light Curve Cleaning and Averaging")
        self.f.write(f"\n\nTimestamp: {timestamp.strftime('%B %d, %Y at %I:%M:%S %p')}")
        full_command = " ".join(sys.argv)
        self.f.write(f"\n\nCommand run: `{full_command}`")

        self.f.write(
            f'\n\nThe SN light curves are separated by filter and labelled as such in the file name. Averaged light curves contain an additional number in the file name that represents the MJD bin size used. Control light curves are located in the "controls" subdirectory and follow the same naming scheme, only with their control index added after the SN name.'
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
        if self.cut_list.has(BadDayCut.name()):
            self.f.write(
                f"\n\t- Averaged light curves (for MJD bin size {format_float_string(mjdbinsize)} days): {self.tnsname}.o.{format_float_string(mjdbinsize)}days.lc.txt and {self.tnsname}.c.{format_float_string(mjdbinsize)}days.lc.txt"
            )
        if self.cut_list.has(ControlLightCurveCut.name()):
            self.f.write(
                f"\n\t- Control light curves, where X=001,...,{num_controls:03d}: {self.tnsname}_iX.o.lc.txt and {self.tnsname}_iX.c.lc.txt"
            )

        self.f.write(
            f'\n\nThe following summarizes the hex values in the "Mask" column of each light curve for each cut applied (see below sections for more information on each cut): '
        )
        if self.cut_list.has(UncertaintyCut.name()):
            self.f.write(
                f"\n\t- Uncertainty cut: {hex(self.cut_list.get(UncertaintyCut.name()).flag)}"
            )
        if self.cut_list.has(ChiSquareCut.name()):
            self.f.write(
                f"\n\t- Chi-square cut: {hex(self.cut_list.get(ChiSquareCut.name()).flag)}"
            )
        if self.cut_list.has(ControlLightCurveCut.name()):
            self.f.write(
                f"\n\t- Control light curve cut: {hex(self.cut_list.get(ControlLightCurveCut.name()).flag)}"
            )
        if badday_cut is not None:
            self.f.write(
                f"\n\t- Bad day (for averaged light curves only): {hex(badday_cut.flag)}"
            )

        custom_cuts = self.cut_list.get_custom_cuts()
        for cut in custom_cuts.values():
            self.f.write(f'\n\t- Custom cut on "{cut.column}" column: {hex(cut.flag)}')

    def add_filter_section(self, filt):
        self.f.write(f"\n\n## FILTER: {filt}")

    def add_template_correction_section(self, output):
        self.f.write(f"\n\n### ATLAS template change correction\n")
        self.f.write("\n".join(output))

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

    def add_uncert_cut_section(self, cut: UncertaintyCut, percent_cut):
        self.f.write(f"\n\n### Uncertainty cut\n")
        self.f.write(f"\nSelected uncertainty max value of {cut.max_value}\n")
        self.f.write(
            f"\nTotal percent of SN light curve flagged with {hex(cut.flag)}: {percent_cut:0.2f}%"
        )

    def add_x2_cut_section(
        self, cut: ChiSquareCut, percent_contamination, percent_loss, percent_cut
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
        cut: ControlLightCurveCut,
        x2_percent_cut,
        stn_percent_cut,
        Nclip_percent_cut,
        Ngood_percent_cut,
        questionable_percent_cut,
        percent_cut,
    ):
        self.f.write(f"\n\n### Control light curve cut\n")
        self.f.write(
            f"\nPercent of SN light curve above x2_max bound ({hex(cut.x2_flag)}): {x2_percent_cut:0.2f}%"
        )
        self.f.write(
            f"\nPercent of SN light curve above snr_max bound ({hex(cut.snr_flag)}): {stn_percent_cut:0.2f}%"
        )
        self.f.write(
            f"\nPercent of SN light curve above Nclip_max bound ({hex(cut.Nclip_flag)}): {Nclip_percent_cut:0.2f}%"
        )
        self.f.write(
            f"\nPercent of SN light curve below Ngood_min bound ({hex(cut.Ngood_flag)}): {Ngood_percent_cut:0.2f}%"
        )
        self.f.write(
            f"\nTotal percent of SN light curve flagged as questionable (not masked with control light curve flags but Nclip > 0) ({hex(cut.questionable_flag)}): {questionable_percent_cut:0.2f}%"
        )
        self.f.write(
            f"\nTotal percent of SN light curve flagged as bad ({hex(cut.flag)}): {percent_cut:0.2f}%"
        )

    def add_badday_cut_section(self, percent_cut):
        self.f.write(
            f'\n\nAfter the cuts are applied, the light curves are resaved with the new "Mask" column.'
        )
        self.f.write(f"\n\n### Bad day cut (averaging)\n")
        self.f.write(
            f"\nTotal percent of SN light curve flagged as bad ({hex(self.cut_list.get_all_default_flags())}): {percent_cut:0.2f}%"
        )
        self.f.write(
            f"\nThe averaged light curves are then saved in a new file with the MJD bin size added to the filename."
        )

    def add_custom_cut_section(self, cut: CustomCut, percent_cut):
        self.f.write(f"\n\n### {cut.name()}\n")
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
            self.t = pd.read_table(self.filename, sep="\s+")
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
                    f"true uncertainties estimation table has {len(matching_ix)} matching rows for TNS name {tnsname} and filter {filt}"
                )

            if len(matching_ix) > 0:
                # update existing row
                idx = matching_ix[0]
                self.t.loc[idx, :] = row
            else:
                self.t = new_row(self.t, row)
        else:
            self.t = new_row(self.t, row)

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
            self.t = pd.read_table(self.filename, sep="\s+")
            print("Success")
        except:
            print(f"No existing chi-square cut table; creating blank table...")
            self.t = pd.DataFrame(
                columns=[
                    "tnsname",
                    "filter",
                    "x2_cut",
                    "use_preSN_lc",
                    "snr_bound",
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
                    f"chi-square cut table has {len(matching_ix)} matching rows for TNS name {tnsname} and filter {filt}"
                )

            if len(matching_ix) > 0:
                # update existing row
                idx = matching_ix[0]
                self.t.loc[idx, :] = row
            else:
                self.t = new_row(self.t, row)
        else:
            self.t = new_row(self.t, row)

    def save(self):
        print(f"\nSaving chi-square cut table at {self.filename}...")
        self.t.to_string(self.filename)
        print("Success")


class CleanLoop:
    def __init__(
        self,
        colnames: PresetColumnNames,
        input_dir: str,
        output_dir: str,
        credentials: Credentials,
        sninfo_filename: Optional[str] = None,
        flux2mag_sigmalimit: float = 3.0,
        overwrite: bool = False,
    ):
        self.colnames = colnames
        self.sn: Optional[Supernova] = None
        self.avg_sn: Optional[AveragedSupernova] = None
        self.cut_list: Optional[CutList] = None
        self.f: Optional[OutputReadMe] = None
        self.p: Optional[PlotPdf] = None

        self.credentials: Credentials = credentials
        self.input_dir: str = input_dir
        self.output_dir: str = output_dir
        self.flux2mag_sigmalimit: float = flux2mag_sigmalimit
        self.overwrite: bool = overwrite

        self.sninfo: SnInfoTable = SnInfoTable(
            self.output_dir, filename=sninfo_filename
        )
        self.uncert_est_info: UncertEstTable = UncertEstTable(self.output_dir)
        if cut_list.has(ChiSquareCut.name()):
            self.x2_cut_info: ChiSquareCutTable = ChiSquareCutTable(self.output_dir)

    def apply_template_correction(
        self,
        maskval=None,
        region1_offset=None,
        region2_offset=None,
        region3_offset=None,
        num_measurements=40,
        plot: bool = False,
    ):
        print(f"\nApplying ATLAS template change correction:")
        if self.sn is None:
            raise RuntimeError("Supernova (self.sn) cannot be None")

        output = self.sn.apply_template_correction(
            maskval=maskval,
            region1_offset=region1_offset,
            region2_offset=region2_offset,
            region3_offset=region3_offset,
            num_measurements=num_measurements,
        )
        print("\n".join(output))

        if self.f is None:
            raise RuntimeError("Output README file (self.f) cannot be None")
        self.f.add_template_correction_section(output)

        if plot:
            if self.p is None:
                raise RuntimeError("Output plots (self.p) cannot be None")
            self.p.plot_template_correction(self.sn.lcs[0])

    def check_uncert_est(
        self, cut: UncertaintyEstimation, apply_function: Callable, plot: bool = False
    ):
        print(f"\nChecking true uncertainties estimation:")
        if self.sn is None:
            raise RuntimeError("Supernova (self.sn) cannot be None")

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
                    if self.p is None:
                        raise RuntimeError("Output plots (self.p) cannot be None")
                    self.p.plot_uncert_est(self.sn)
        else:
            print("True uncertainties estimation not needed; skipping procedure...")

        if self.f is None:
            raise RuntimeError("Output README file (self.f) cannot be None")
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

    def apply_uncert_cut(self, cut: UncertaintyCut | None, plot: bool = False):
        if cut is None:
            return

        print(f"\nApplying uncertainty cut ({cut}):")
        if self.sn is None:
            raise RuntimeError("Supernova (self.sn) cannot be None")
        percent_cut = self.sn.apply_cut(cut)
        print("Success")
        print(
            f"Total percent of SN light curve flagged with {hex(cut.flag)}: {percent_cut:0.2f}%"
        )

        if self.f is None:
            raise RuntimeError("Output README file (self.f) cannot be None")
        self.f.add_uncert_cut_section(cut, percent_cut)

        if plot:
            if self.p is None:
                raise RuntimeError("Output plots (self.p) cannot be None")
            self.p.plot_cut(self.sn, cut.flag, title=UncertaintyCut.name())

    def apply_x2_cut(self, cut: ChiSquareCut | None, plot: bool = False):
        if cut is None:
            return None

        print(f"\nApplying chi-square cut ({cut}):")
        if self.sn is None:
            raise RuntimeError("Supernova (self.sn) cannot be None")
        if self.sn.colnames.chisquare is None:
            print(
                "WARNING: No chi-square column name provided in config file; skipping..."
            )
            return None

        if cut.use_pre_mjd0_lc:
            print("Using pre-MJD0 light curve to determine contamination and loss")
            if self.sn.mjd0 is None:
                raise RuntimeError(
                    "MJD0 cannot be None. Please provide MJD0 throught the SN info table or --mjd0 argument, or set the use_pre_MJD0_lc field in the config file to False."
                )
            lc_temp = deepcopy(self.sn.lcs[0])
            ix = lc_temp.ix_inrange(lc_temp.colnames.mjd, uplim=self.sn.mjd0)
            if len(ix) < 1:
                raise RuntimeError("no pre-MJD0 light curve available")
        else:
            print("Using control light curves to determine contamination and loss")
            if self.sn.num_controls < 1:
                raise RuntimeError(
                    "No control light curves loaded. Use the --num_controls argument to load control light curves, or change the [x2_cut][use_pre_mjd0_lc] field to True."
                )
            lc_temp = self.sn.get_all_controls()
            ix = lc_temp.t.index.values

        limcuts = LimCutsTable(lc_temp, cut.snr_bound, indices=ix)
        limcuts.calculate_table(cut.min_cut, cut.max_cut, cut.cut_step)
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
            if self.p is None:
                raise RuntimeError("Output plots (self.p) cannot be None")
            self.p.plot_limcuts(limcuts, cut)
            self.p.plot_cut(self.sn, cut.flag, title=ChiSquareCut.name())

        if self.f is None:
            raise RuntimeError("Output README file (self.f) cannot be None")
        self.f.add_x2_cut_section(
            cut, data["Pcontamination"], data["Ploss"], percent_cut
        )

        x2_info_row = {
            "tnsname": self.sn.tnsname,
            "filter": self.sn.filt,
            "x2_cut": cut.max_value,
            "use_pre_mjd0_lc": cut.use_pre_mjd0_lc,
            "snr_bound": cut.snr_bound,
            "pct_contamination": round(data["Pcontamination"], 2),
            "pct_loss": round(data["Ploss"], 2),
        }
        return x2_info_row

    def apply_controls_cut(
        self, cut: ControlLightCurveCut | None, previous_flags: int, plot: bool = False
    ):
        if cut is None:
            return

        print(f"\nApplying control light curve cut ({cut}):")
        if self.sn is None:
            raise RuntimeError("Supernova (self.sn) cannot be None")

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
            f"Percent of data above x2_max bound ({hex(cut.x2_flag)}): {x2_percent_cut:0.2f}%"
        )
        print(
            f"Percent of data above snr_max bound ({hex(cut.snr_flag)}): {stn_percent_cut:0.2f}%"
        )
        print(
            f"Percent of data above Nclip_max bound ({hex(cut.Nclip_flag)}): {Nclip_percent_cut:0.2f}%"
        )
        print(
            f"Percent of data below Ngood_min bound ({hex(cut.Ngood_flag)}): {Ngood_percent_cut:0.2f}%"
        )
        print(
            f"Total percent of data flagged as questionable (not masked with control light curve flags but Nclip > 0) ({hex(cut.questionable_flag)}): {questionable_percent_cut:0.2f}%"
        )
        print(
            f"Total percent of data flagged as bad ({hex(cut.flag)}): {percent_cut:0.2f}%"
        )

        if self.f is None:
            raise RuntimeError("Output README file (self.f) cannot be None")
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
            if self.p is None:
                raise RuntimeError("Output plots (self.p) cannot be None")
            self.p.plot_cut(self.sn, cut.flag, title=ControlLightCurveCut.name())

    def apply_badday_cut(
        self, cut: BadDayCut | None, previous_flags, plot: bool = False
    ):
        if cut is None:
            return

        print(
            f"\nApplying bad day cut (averaging) with MJD bin size of {cut.mjd_bin_size} days ({cut}):"
        )
        if self.sn is None:
            raise RuntimeError("Supernova (self.sn) cannot be None")

        self.avg_sn, percent_cut = self.sn.apply_badday_cut(
            cut, previous_flags, flux2mag_sigmalimit=self.flux2mag_sigmalimit
        )
        print("Success")
        print(
            f"Total percent of SN light curve flagged as bad ({hex(self.cut_list.get_all_default_flags())}): {percent_cut:0.2f}"
        )

        if self.f is None:
            raise RuntimeError("Output README file (self.f) cannot be None")
        self.f.add_badday_cut_section(percent_cut)

        if plot:
            if self.p is None:
                raise RuntimeError("Output plots (self.p) cannot be None")
            self.p.plot_cut(self.avg_sn, cut.flag, title=BadDayCut.name())
            self.p.plot_averaged_SN(
                self.avg_sn, cut.flag, plot_controls=True, plot_flagged=False
            )

    def apply_custom_cut(self, cut: CustomCut, plot: bool = False):
        print(f"\nApplying custom cut ({cut})...")
        if self.sn is None:
            raise RuntimeError("Supernova (self.sn) cannot be None")

        percent_cut = self.sn.apply_cut(cut)
        print("Success")
        print(
            f"Total percent of SN light curve flagged with {hex(cut.flag)}: {percent_cut:0.2f}%"
        )

        if self.f is None:
            raise RuntimeError("Output README file (self.f) cannot be None")
        self.f.add_custom_cut_section(cut, percent_cut)

        if plot:
            if self.p is None:
                raise RuntimeError("Output plots (self.p) cannot be None")
            self.p.plot_cut(self.sn, cut.flag, title=cut.name())

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
        self.sn = Supernova(self.colnames, tnsname=tnsname, mjd0=mjd0, filt=filt)
        try:
            self.sn.load_all(self.input_dir, num_controls=num_controls)
        except Exception as e:
            raise RuntimeError(f"Could not load light curves: {str(e)}")

        # prepare the light curves for cleaning
        print()
        self.sn.prep_for_cleaning(verbose=True)

        if plot:
            # initialize PDF of diagnostic plots
            self.p = PlotPdf(f"{self.output_dir}/{tnsname}", tnsname, filt=filt)

            # plot original SN light curve and control light curves
            self.p.plot_SN(self.sn, plot_controls=True, plot_template_changes=True)

        if self.cut_list is None:
            raise RuntimeError("CutList (self.cut_list) cannot be None")

        # template correction
        if apply_template_correction:
            self.apply_template_correction()

        # uncertainty cut
        self.apply_uncert_cut(self.cut_list.get(UncertaintyCut.name()), plot=plot)

        # true uncertainties estimation
        _, uncert_est_info_row = self.check_uncert_est(
            self.cut_list.get(UncertaintyEstimation.name()),
            apply_function=apply_uncert_est_function,
            plot=plot,
        )
        self.uncert_est_info.add_row(uncert_est_info_row)

        # chi-square cut
        x2_info_row = self.apply_x2_cut(
            self.cut_list.get(ChiSquareCut.name()), plot=plot
        )
        if x2_info_row is not None and self.cut_list.has(ChiSquareCut.name()):
            self.x2_cut_info.add_row(x2_info_row)

        # control light curve cut
        self.apply_controls_cut(
            self.cut_list.get(ControlLightCurveCut.name()),
            previous_flags=self.cut_list.get_previous_flags(
                ControlLightCurveCut.name()
            ),
            plot=plot,
        )

        # custom cuts
        custom_cuts = self.cut_list.get_custom_cuts()
        for cut in custom_cuts.values():
            self.apply_custom_cut(cut, plot=plot)

        # plot the cleaned light curves so far
        previous_flags = self.cut_list.get_previous_flags(BadDayCut.name())
        if plot and previous_flags > 0:
            print()
            self.p.plot_cut(self.sn, previous_flags, title="All previous cuts")
            self.p.plot_cleaned_SN(
                self.sn, previous_flags, plot_controls=True, plot_flagged=False
            )

        # bad day cut (averaging)
        self.apply_badday_cut(
            self.cut_list.get(BadDayCut.name()),
            previous_flags=previous_flags,
            plot=plot,
        )

        # save cleaned SN and control light curves
        self.sn.save_all(self.output_dir, overwrite=self.overwrite)

        if self.cut_list.has(BadDayCut.name()):
            if self.avg_sn is None:
                raise RuntimeError("Averaged supernova (self.avg_sn) cannot be None")
            # save averaged SN and control light curves
            self.avg_sn.save_all(self.output_dir, overwrite=self.overwrite)

        if self.cut_list.has(UncertaintyEstimation.name()):
            # save uncertainty estimation table
            self.uncert_est_info.save()

        if self.cut_list.has(ChiSquareCut.name()):
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
        filts: Optional[List[str]] = None,
        cut_list: Optional[CutList] = None,
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
                plot
                or (
                    cut_list.has(ChiSquareCut.name())
                    and cut_list.get(ChiSquareCut.name()).use_pre_mjd0_lc
                )
            ):
                mjd0, coords = get_mjd0_from_tns(tnsname, self.sninfo, self.credentials)
                if not coords is None:
                    print(f"Setting MJD0 to {mjd0}")
                    self.sninfo.update_row(tnsname, coords=coords, mjd0=mjd0)
                    print("Success")
            else:
                print(f"\nSetting MJD0 to {mjd0}")

            if filts is None:
                print("Searching for filters in input directory...")
                filts = find_all_filts(self.input_dir, tnsname)
                print(f"Filters found: {filts}")

            for filt in filts:
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


def parse_config_cuts(args, config, colnames):
    cut_list = CutList()
    if args.custom_cuts:
        config_custom_cuts = get_config_custom_cuts(config)

    print(f"\nProcedures to apply:")

    # always check true uncertainties estimation, but will only apply if args.true_uncert_est
    temp_x2_max_value = float(config["uncert_est"]["temp_x2_max_value"])
    print(
        f"- True uncertainties estimation check (using temporary chi-square cut at {temp_x2_max_value})"
    )
    uncert_est = UncertaintyEstimation(
        temp_x2_max_value,
        hexstring_to_int(config["uncert_cut"]["flag"]),
    )
    cut_list.add(uncert_est)

    if args.uncert_cut:
        uncert_cut = UncertaintyCut(
            colnames.dflux,
            flag=hexstring_to_int(config["uncert_cut"]["flag"]),
            max_value=float(config["uncert_cut"]["max_value"]),
        )
        cut_list.add(uncert_cut)
        print(f"- {uncert_cut}")

    if args.x2_cut:
        x2_cut = ChiSquareCut(
            colnames.chisquare,
            flag=hexstring_to_int(config["x2_cut"]["flag"]),
            max_value=float(config["x2_cut"]["max_value"]),
            snr_bound=float(config["x2_cut"]["snr_bound"]),
            min_cut=int(config["x2_cut"]["min_cut"]),
            max_cut=int(config["x2_cut"]["max_cut"]),
            cut_step=int(config["x2_cut"]["cut_step"]),
            use_pre_mjd0_lc=parse_config_str(config["x2_cut"]["use_pre_mjd0_lc"]),
        )
        cut_list.add(x2_cut)
        print(f"- {x2_cut}")

    if args.controls_cut:
        controls_cut = ControlLightCurveCut(
            flag=hexstring_to_int(config["controls_cut"]["bad_flag"]),
            questionable_flag=hexstring_to_int(
                config["controls_cut"]["questionable_flag"]
            ),
            x2_max=float(config["controls_cut"]["x2_max"]),
            x2_flag=hexstring_to_int(config["controls_cut"]["x2_flag"]),
            snr_max=float(config["controls_cut"]["snr_max"]),
            snr_flag=hexstring_to_int(config["controls_cut"]["snr_flag"]),
            Nclip_max=int(config["controls_cut"]["Nclip_max"]),
            Nclip_flag=hexstring_to_int(config["controls_cut"]["Nclip_flag"]),
            Ngood_min=int(config["controls_cut"]["Ngood_min"]),
            Ngood_flag=hexstring_to_int(config["controls_cut"]["Ngood_flag"]),
        )
        cut_list.add(controls_cut)
        print(f"- {controls_cut}")

    if args.averaging:
        badday_cut = BadDayCut(
            flag=hexstring_to_int(config["averaging"]["flag"]),
            mjd_bin_size=(
                float(config["averaging"]["mjd_bin_size"])
                if args.mjd_bin_size is None
                else args.mjd_bin_size
            ),
            x2_max=float(config["averaging"]["x2_max"]),
            Nclip_max=int(config["averaging"]["Nclip_max"]),
            Ngood_min=int(config["averaging"]["Ngood_min"]),
            ixclip_flag=hexstring_to_int(config["averaging"]["ixclip_flag"]),
            smallnum_flag=hexstring_to_int(config["averaging"]["smallnum_flag"]),
        )
        cut_list.add(badday_cut)
        print(f"- Averaging / {badday_cut}")

    if args.custom_cuts:
        for i in range(len(config_custom_cuts)):
            cut_settings = config_custom_cuts[i]
            try:
                custom_cut = CustomCut(
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
                cut_list.add(custom_cut)
                print(f"- Custom cut {i}: {custom_cut.name()}")
            except Exception as e:
                print(
                    f"WARNING: Could not parse custom cut {i}: {cut_settings}. Error: {str(e)}"
                )

    duplicate_flags = cut_list.get_flag_duplicates()
    if len(duplicate_flags) > 0:
        raise ValueError(
            f"Cuts in the config file contain duplicate flags: {duplicate_flags}."
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
        "-p",
        "--preset",
        type=str,
        default="atlas",
        help="preset name from config file (ex. atlas, rubin, tess)",
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
        raise RuntimeError("Please specify at least one TNS name to clean.")
    if len(args.tnsnames) > 1 and not args.mjd0 is None:
        raise RuntimeError(f"Cannot specify one MJD0 {args.mjd0} for a batch of SNe.")
    print(f"\nList of transients to clean: {args.tnsnames}")

    allowed_presets = get_allowed_presets(config)
    if args.preset is None or args.preset not in allowed_presets:
        raise RuntimeError(
            f"Please specify the preset name to load from the config file (allowed presets: {allowed_presets})"
        )
    print(f"\nLoading '{args.preset}' preset column names from config.ini...")
    colnames = PresetColumnNames(config, args.preset)
    print(colnames.__str__())
    print("Success")

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
    filters = parse_comma_separated_string(args.filters)
    if filters is not None:
        print(f"Filters to clean: {filters}")
    flux2mag_sigmalimit = float(config["download"]["flux2mag_sigmalimit"])
    print(f"Sigma limit when converting flux to magnitude: {flux2mag_sigmalimit}")
    if args.mjd0:
        print(f"MJD0: {args.mjd0}")

    # print(f'\nApplyin control light curve cut: {args.controls}')
    num_controls = (
        args.num_controls
        if not args.num_controls is None
        else int(config["download"]["num_controls"])
    )
    print(f"Number of control light curves to clean: {num_controls}")

    cut_list = parse_config_cuts(args, config, colnames)

    print()
    credentials = Credentials(
        config["credentials"]["atlas_username"],
        config["credentials"]["atlas_password"],
        config["credentials"]["tns_api_key"],
        config["credentials"]["tns_id"],
        config["credentials"]["tns_bot_name"],
    )
    clean = CleanLoop(
        colnames,
        input_dir,
        output_dir,
        credentials,
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
        filts=filters,
        apply_template_correction=args.template_correction,
        plot=args.plot,
    )
