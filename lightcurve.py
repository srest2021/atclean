#!/usr/bin/env python

from abc import ABC, abstractmethod
from configparser import ConfigParser
from typing import Dict, Any, List, Optional, Self, Set, Tuple, Type
import re, json, requests, time, sys, io, bisect
from astropy import units as u
from astropy.coordinates import Angle
from astropy.time import Time
from collections import OrderedDict

import scipy
from pdastro import AnotB, pdastrostatsclass
import numpy as np
import pandas as pd
from copy import deepcopy
from pathlib import Path
from utils import (
    DISC_DATE_BUFFER,
    TEMPLATE_CHANGE_1_MJD,
    TEMPLATE_CHANGE_2_MJD,
    AandB,
    AorB,
    BadDayCut,
    ControlLightCurveCut,
    Coordinates,
    Cut,
    PresetColumnNames,
    UncertaintyEstimation,
    combine_flags,
    find_all_control_indices,
    format_float,
    get_filename,
    get_tns_coords_from_json,
    get_tns_mjd0_from_json,
    new_row,
    query_atlas,
    query_tns,
    PlotLimits,
)


class Supernova:
    def __init__(
        self,
        colnames: PresetColumnNames,
        tnsname: str = None,
        ra: str = None,
        dec: str = None,
        mjd0: float = None,
        filt: str = "o",
    ):
        self.colnames = deepcopy(colnames)

        self.tnsname = tnsname
        self.coords: Coordinates = Coordinates(ra, dec)
        self.mjd0 = mjd0
        self.filt = filt

        self.lcs: Dict[int, LightCurve] = {}

        self.num_controls = 0
        self._all_indices = None
        self._control_indices = None

    def get(self, control_index=0):
        try:
            return self.lcs[control_index].t
        except:
            raise RuntimeError(
                f"Cannot get control light curve {control_index}. Num controls set to {self.num_controls} and {len(self.lcs)} lcs in dictionary."
            )

    def get_lims(self, control_index=0, flag: Optional[int] = None):
        lims = PlotLimits()

        lims.set_xlims(
            self.lcs[control_index].get_xlims(
                mjd0=self.mjd0 if control_index == 0 else None
            )
        )

        lims.set_ylims(self.lcs[control_index].get_ylims(flag=flag))

        return lims

    def get_tns_data(self, api_key, tns_id, bot_name):
        if self.coords.is_incomplete() or self.mjd0 is None:
            print(f"\nQuerying TNS for {self.tnsname} data...")
            json_data = query_tns(self.tnsname, api_key, tns_id, bot_name)
            if json_data is None:
                print(f"Skipping...")
                return

            if self.coords.is_incomplete():
                self.coords = get_tns_coords_from_json(json_data)

            if self.mjd0 is None:
                self.mjd0 = get_tns_mjd0_from_json(json_data)

            print("Success")

    def verify_mjds(self, verbose=False):
        """Sort SN and control light curves by MJD"""
        self.lcs[0].t.sort_values(
            by=[self.colnames.mjd], ignore_index=True, inplace=True
        )

        if self.num_controls == 0:
            return

        if verbose:
            print("\nMaking sure SN and control light curve MJDs match up exactly:")

        sn_sorted_mjd = self.lcs[0].t[self.colnames.mjd].to_numpy()

        for control_index in self.control_lc_indices:
            # sort by MJD
            self.lcs[control_index].t.sort_values(
                by=[self.colnames.mjd], ignore_index=True, inplace=True
            )
            control_sorted_mjd = self.lcs[control_index].t[self.colnames.mjd].to_numpy()

            if (len(sn_sorted_mjd) != len(control_sorted_mjd)) or not np.array_equal(
                sn_sorted_mjd, control_sorted_mjd
            ):
                if verbose:
                    print(
                        f"MJDs out of agreement for control light curve {control_index}, fixing..."
                    )

                only_sn_mjd = AnotB(sn_sorted_mjd, control_sorted_mjd)
                only_control_mjd = AnotB(control_sorted_mjd, sn_sorted_mjd)

                # for the MJDs only in SN, add row with that MJD to control light curve,
                # with all values of other columns NaN
                if len(only_sn_mjd) > 0:
                    for mjd in only_sn_mjd:
                        self.lcs[control_index].newrow(
                            {
                                self.colnames.mjd: mjd,
                                self.colnames.mask: 0,
                            }
                        )

                # remove indices of rows in control light curve for which there is no MJD in the SN lc
                if len(only_control_mjd) > 0:
                    ix_to_skip = []
                    for mjd in only_control_mjd:
                        matching_ix = self.lcs[control_index].ix_equal(
                            self.colnames.mjd, mjd
                        )
                        if len(matching_ix) != 1:
                            raise RuntimeError(
                                f"Couldn't find MJD={mjd} in MJD column, but should be there!"
                            )
                        ix_to_skip.extend(matching_ix)
                    ix = AnotB(self.lcs[control_index].getindices(), ix_to_skip)
                else:
                    ix = self.lcs[control_index].getindices()

                # sort again
                sorted_ix = self.lcs[control_index].ix_sort_by_cols(
                    self.colnames.mjd, indices=ix
                )
                self.lcs[control_index].t = self.lcs[control_index].t.loc[sorted_ix]

            self.lcs[control_index].t.reset_index(drop=True, inplace=True)

        print("Success")

    def prep_for_cleaning(self, verbose=False):
        if verbose:
            print(
                f"Adding blank '{self.colnames.mask}' columns, replacing infs with NaNs, and calculating flux/dflux..."
            )

        for control_index in self.lc_indices:
            # add blank 'Mask' column
            self.lcs[control_index].t[self.colnames.mask] = 0
            # remove rows with duJy=0 or uJy=NaN
            self.lcs[control_index].remove_invalid_rows()
            # calculate flux/dflux column
            self.lcs[control_index].calculate_fdf_column()
        print("Success")

        # make sure SN and control lc MJDs match up exactly
        self.verify_mjds(verbose=verbose)

    def apply_template_correction(
        self,
        maskval=None,
        region1_offset=None,
        region2_offset=None,
        region3_offset=None,
        num_measurements=40,
    ):
        if self.mjd0 is None:
            raise RuntimeError("Cannot apply template correction without MJD0")
        return self.lcs[0].apply_template_correction(
            self.mjd0,
            maskval=maskval,
            region1_offset=region1_offset,
            region2_offset=region2_offset,
            region3_offset=region3_offset,
            num_measurements=num_measurements,
        )

    def apply_cut(self, cut: Cut):
        sn_percent_cut = None
        for control_index in self.lc_indices:
            percent_cut = self.lcs[control_index].apply_cut(cut)
            if control_index == 0:
                sn_percent_cut = percent_cut
        return sn_percent_cut

    def get_uncert_est_stats(self, cut: UncertaintyEstimation):
        def get_sigma_extra(median_dflux, stdev):
            diff = stdev**2 - median_dflux**2
            return max(0, np.sqrt(diff)) if diff > 0 else 0

        stats = pd.DataFrame(
            columns=["control_index", "median_dflux", "stdev", "sigma_extra"]
        )
        stats["control_index"] = self.control_lc_indices
        stats.set_index("control_index", inplace=True)

        use_x2_clean_ix = (
            self.colnames.chisquare is not None
            and self.colnames.chisquare in self.lcs[0].t.columns
        )

        for control_index in self.control_lc_indices:
            dflux_clean_ix = self.lcs[control_index].ix_unmasked(
                self.colnames.mask, maskval=cut.uncert_cut_flag
            )

            if use_x2_clean_ix:
                x2_clean_ix = self.lcs[control_index].ix_inrange(
                    colnames=[self.colnames.chisquare],
                    uplim=cut.temp_x2_max_value,
                    exclude_uplim=True,
                )
                clean_ix = AandB(dflux_clean_ix, x2_clean_ix)
            else:
                clean_ix = dflux_clean_ix

            median_dflux = self.lcs[control_index].get_median_dflux(indices=clean_ix)

            stdev_flux = self.lcs[control_index].get_stdev_flux(indices=clean_ix)
            if stdev_flux is None:
                print(
                    f'WARNING: Could not get flux std dev using clean indices; retrying without preliminary chi-square cut of {cut.params["temp_x2_max_value"]}...'
                )
                stdev_flux = self.lcs[control_index].get_stdev_flux(
                    indices=dflux_clean_ix
                )
                if stdev_flux is None:
                    print(
                        "WARNING: Could not get flux std dev using clean indices; retrying with all indices..."
                    )
                    stdev_flux = self.lcs[control_index].get_stdev_flux(
                        control_index=control_index
                    )

            sigma_extra = get_sigma_extra(median_dflux, stdev_flux)

            stats.loc[control_index, "median_dflux"] = median_dflux
            stats.loc[control_index, "stdev"] = stdev_flux
            stats.loc[control_index, "sigma_extra"] = sigma_extra

        return stats

    def add_noise_to_dflux(self, sigma_extra):
        for control_index in self.lc_indices:
            self.lcs[control_index].add_noise_to_dflux(sigma_extra)

    def get_all_controls(self):
        controls = [
            deepcopy(self.lcs[control_index].t)
            for control_index in self.lcs
            if control_index > 0
        ]
        all_controls = LightCurve(self.colnames)
        all_controls.t = pd.concat(controls, ignore_index=True)
        return all_controls

    def calculate_control_stats(self, previous_flags):
        print("Calculating control light curve statistics...")

        len_mjd = len(self.lcs[0].t[self.colnames.mjd])

        # construct arrays for control lc data
        uJy = np.full((self.num_controls, len_mjd), np.nan)
        duJy = np.full((self.num_controls, len_mjd), np.nan)
        Mask = np.full((self.num_controls, len_mjd), 0, dtype=np.int32)

        i = 1
        for control_index in self.control_lc_indices:
            if len(self.lcs[control_index].t) != len_mjd or not np.array_equal(
                self.lcs[0].t[self.colnames.mjd],
                self.lcs[control_index].t[self.colnames.mjd],
            ):
                raise RuntimeError(
                    f"SN lc not equal to control lc for control_index {control_index}! Rerun or debug verify_mjds()."
                )
            else:
                uJy[i - 1, :] = self.lcs[control_index].t[self.colnames.flux]
                duJy[i - 1, :] = self.lcs[control_index].t[
                    self.lcs[control_index].colnames.dflux_new
                ]
                Mask[i - 1, :] = self.lcs[control_index].t[self.colnames.mask]

            i += 1

        c2_param2columnmapping = self.lcs[0].intializecols4statparams(
            prefix="c2_", format4outvals="{:.2f}", skipparams=["converged", "i"]
        )

        for index in range(uJy.shape[-1]):
            pda4MJD = pdastrostatsclass()
            pda4MJD.t[self.colnames.flux] = uJy[0:, index]
            pda4MJD.t[self.lcs[0].colnames.dflux_new] = duJy[0:, index]
            pda4MJD.t[self.colnames.mask] = np.bitwise_and(
                Mask[0:, index], previous_flags
            )

            pda4MJD.calcaverage_sigmacutloop(
                self.colnames.flux,
                noisecol=self.lcs[0].colnames.dflux_new,
                maskcol=self.colnames.mask,
                maskval=previous_flags,
                verbose=1,
                Nsigma=3.0,
                median_firstiteration=True,
            )
            self.lcs[0].statresults2table(
                pda4MJD.statparams, c2_param2columnmapping, destindex=index
            )

    def apply_controls_cut(self, cut: ControlLightCurveCut, previous_flags: int):
        self.calculate_control_stats(previous_flags)
        self.lcs[0].t["c2_abs_stn"] = (
            self.lcs[0].t["c2_mean"] / self.lcs[0].t["c2_mean_err"]
        )

        # flag SN measurements
        self.lcs[0].flag_by_control_stats(cut)

        # copy over SN's control cut flags to control light curve 'Mask' columns
        flags_arr = np.full(
            self.lcs[0].t[self.colnames.mask].shape,
            combine_flags(cut.get_flags()),
        )
        flags_to_copy = np.bitwise_and(self.lcs[0].t[self.colnames.mask], flags_arr)
        for control_index in self.control_lc_indices:
            self.lcs[control_index].copy_flags(flags_to_copy)

        # self.drop_extra_columns()

        len_ix = len(self.lcs[0].getindices())
        x2_percent_cut = (
            100
            * len(self.lcs[0].ix_masked(self.colnames.mask, maskval=cut.x2_flag))
            / len_ix
        )
        stn_percent_cut = (
            100
            * len(self.lcs[0].ix_masked(self.colnames.mask, maskval=cut.snr_flag))
            / len_ix
        )
        Nclip_percent_cut = (
            100
            * len(self.lcs[0].ix_masked(self.colnames.mask, maskval=cut.Nclip_flag))
            / len_ix
        )
        Ngood_percent_cut = (
            100
            * len(self.lcs[0].ix_masked(self.colnames.mask, maskval=cut.Ngood_flag))
            / len_ix
        )
        questionable_percent_cut = (
            100
            * len(
                self.lcs[0].ix_masked(self.colnames.mask, maskval=cut.questionable_flag)
            )
            / len_ix
        )
        percent_cut = (
            100
            * len(self.lcs[0].ix_masked(self.colnames.mask, maskval=cut.flag))
            / len_ix
        )
        return (
            x2_percent_cut,
            stn_percent_cut,
            Nclip_percent_cut,
            Ngood_percent_cut,
            questionable_percent_cut,
            percent_cut,
        )

    def apply_badday_cut(self, cut: BadDayCut, previous_flags, flux2mag_sigmalimit=3.0):
        avg_sn = AveragedSupernova(
            self.colnames,
            tnsname=self.tnsname,
            mjd0=self.mjd0,
            filt=self.filt,
            mjdbinsize=cut.mjd_bin_size,
            flag=cut.flag,
        )
        avg_sn.num_controls = self.num_controls

        for control_index in self.lc_indices:
            avg_sn.set_avg_lc(
                self.lcs[control_index].average(
                    cut,
                    previous_flags,
                    mjdbinsize=cut.mjd_bin_size,
                    flux2mag_sigmalimit=flux2mag_sigmalimit,
                ),
                control_index=control_index,
            )

        all_flags = previous_flags | combine_flags(cut.get_flags())
        percent_cut = (
            100
            * len(avg_sn.lcs[0].ix_masked(self.colnames.mask, maskval=all_flags))
            / len(avg_sn.lcs[0].t)
        )
        return avg_sn, percent_cut

    def drop_extra_columns(self):
        for control_index in self.lc_indices:
            self.lcs[control_index].drop_extra_columns()

    def count_files_in_dir(self, path):
        directory_path = Path(path)
        files = [f for f in directory_path.iterdir() if f.is_file()]
        return len(files)

    def remove_flag(self, flag):
        for control_index in self.lc_indices:
            self.lcs[control_index].remove_flag(flag)

    def load(self, input_dir, control_index=0, cleaned=False):
        self.lcs[control_index] = LightCurve(
            self.colnames, control_index=control_index, filt=self.filt
        )
        self.lcs[control_index].load_lc(input_dir, self.tnsname, cleaned=cleaned)

    def load_all(self, input_dir, num_controls=0, cleaned=False):
        self.lcs = {}
        self.num_controls = 0

        print(f"\nLoading SN light curve and {num_controls} control light curves...")

        # load SN light curve
        self.load(input_dir, cleaned=cleaned)

        control_indices = find_all_control_indices(
            input_dir, self.tnsname, filt=self.filt
        )
        if len(control_indices) < num_controls:
            raise RuntimeError(
                f"Tried to load {num_controls} control light curves, but only {len(control_indices)} found: {control_indices}"
            )

        if num_controls > 0:
            # keep iterating over control indices until we successfully load num_controls light curves
            control_index = 1
            while self.num_controls < num_controls:
                try:
                    self.load(input_dir, control_index=control_index, cleaned=cleaned)
                    self.num_controls += 1
                except:
                    print(
                        f"Could not load control light curve {control_index}; skipping..."
                    )
                    del self.lcs[control_index]
                control_index += 1

        print(
            f"Successfully loaded SN light curve and {self.num_controls} control light curves (control indices: {self.control_lc_indices})"
        )

        # check for dflux_new column if cleaned
        # if found, update colnames in self and all lc objects
        if cleaned and f"{self.colnames.dflux}_new" in self.lcs[0].t.columns:
            self.update_all_colnames("dflux_new", f"{self.colnames.dflux}_new")

    def update_all_colnames(self, key: str, name: str):
        print(
            f"Updating column names for all light curves in this Supernova object (key: {key}, name: {name})..."
        )
        self.colnames.update(key, name)
        for control_index in self.lc_indices:
            self.lcs[control_index].colnames.update(key, name)

    @property
    def lc_indices(self):
        if not self._all_indices:
            self._all_indices = list(self.lcs.keys())
            self._all_indices.sort()
        return self._all_indices

    @property
    def control_lc_indices(self):
        if not self._control_indices:
            self._control_indices = list(self.lcs.keys())
            if 0 in self._control_indices:
                self._control_indices.remove(0)
            self._control_indices.sort()
        return self._control_indices

    def remove_lc_index(self, index: int):
        if index not in self.lcs.keys():
            raise ValueError(
                f"Cannot remove control index {index} because there is no such light curve"
            )
        if index not in self.lc_indices or index not in self.control_lc_indices:
            print(
                f"WARNING: Cannot remove control index {index} because it has already been removed"
            )
            return
        if index == 0:
            raise ValueError(
                "Cannot remove control index 0 because it is reserved for the SN light curve"
            )

        self._all_indices.remove(index)
        self._control_indices.remove(index)

    def add_lc_index(self, index: int):
        if index not in self.lcs.keys():
            raise ValueError(
                f"Cannot add control index {index} because there is no such light curve"
            )
        if index in self.lc_indices or index in self.control_lc_indices:
            print(
                f"WARNING: Cannot add control index {index} because it has already been added"
            )
            return

        bisect.insort(self._all_indices, index)
        bisect.insort(self._control_indices, index)

    def remove_lc_indices(self, indices: List[int]):
        for index in indices:
            self.remove_lc_index(index)

    def add_lc_indices(self, indices: List[int]):
        for index in indices:
            self.add_lc_index(index)

    def reset_lc_indices(self):
        for index in self.lcs.keys():
            if index not in self.lc_indices:
                bisect.insort(self._all_indices, index)
            if index not in self.control_lc_indices:
                bisect.insort(self._control_indices, index)

    def save_all(self, output_dir, overwrite=False, cleaned=True):
        print(
            f'\nDropping extra columns and saving {"cleaned " if cleaned else ""}SN light curve and {self.num_controls} {"cleaned " if cleaned else ""}control light curves...'
        )
        for control_index in self.lc_indices:
            self.lcs[control_index].drop_extra_columns()
            self.lcs[control_index].save_lc(
                output_dir, self.tnsname, overwrite=overwrite, cleaned=cleaned
            )
        print("Success")

    def __str__(self):
        return f"SN {self.tnsname} at {self.coords}: MJD0 = {self.mjd0}, {self.num_controls} control light curves"


class AveragedSupernova(Supernova):
    def __init__(
        self,
        colnames: PresetColumnNames,
        tnsname: str = None,
        ra: str = None,
        dec: str = None,
        mjd0: float | None = None,
        mjdbinsize: float = 1.0,
        filt: str = "o",
        flag: int = 0x800000,
    ):
        Supernova.__init__(self, colnames, tnsname, ra, dec, mjd0, filt)
        self.flag = flag
        self.mjdbinsize = mjdbinsize
        self._mjd_ranges = None

        self.lcs: Dict[int, AveragedLightCurve] = {}

    def set_avg_lc(self, lc, control_index: int = 0):
        self.lcs[control_index] = deepcopy(lc)

    def set_avg_lcs(self, lcs: Dict):
        self.lcs = deepcopy(lcs)

    def get_avg(self, control_index: int = 0):
        try:
            return self.lcs[control_index].t
        except:
            raise RuntimeError(
                f"Cannot get averaged control light curve {control_index}. Num controls set to {self.num_controls} and {len(self.lcs)} lcs in dictionary."
            )

    def has_pre_mjd0_ix(self):
        for control_index in self.control_lc_indices:
            if not self.lcs[control_index].has_pre_mjd0_ix():
                return False
        return True

    def has_valid_mjd_ix(self):
        for control_index in self.control_lc_indices:
            if not self.lcs[control_index].has_valid_mjd_ix():
                return False
        return True

    def get_good_indices(self, control_index: int = 0, flag: Optional[int] = None):
        # if flag is 0 or no mask column, return all indices
        if flag == 0 or not self.colnames.mask in self.lcs[control_index].t.columns:
            return self.lcs[control_index].getindices()

        if flag is None:
            flag = self.flag

        if flag is None:
            # get all flags present in mask column
            flag = self.lcs[control_index].get_flags()

        return self.lcs[control_index].ix_unmasked(self.colnames.mask, maskval=flag)

    def get_bad_indices(self, control_index: int = 0, flag: Optional[int] = None):
        # if flag is 0 or no mask column, return no indices
        if flag == 0 or not self.colnames.mask in self.lcs[control_index].t.columns:
            return []

        if flag is None:
            flag = self.flag

        if flag is None:
            # get all flags present in mask column
            flag = self.lcs[control_index].get_flags()

        return self.lcs[control_index].ix_masked(self.colnames.mask, maskval=flag)

    def set_mjd_ranges(self, mjd_ranges: List[List[float]]):
        if self._mjd_ranges is not None and set(map(tuple, mjd_ranges)) == set(
            map(tuple, self._mjd_ranges)
        ):
            return

        mjd_ranges.sort()
        self._mjd_ranges = mjd_ranges
        for control_index in self.lc_indices:
            self.lcs[control_index].set_valid_mjd_ix(self._mjd_ranges)

    def set_pre_MJD0_ix(self, mjd0: Optional[float] = None, control_index: int = 0):
        if self.mjd0 is None and mjd0 is None:
            raise RuntimeError("Cannot set pre-MJD0 indices without MJD0")

        if mjd0 is not None:
            self.mjd0 = mjd0

        self.lcs[control_index].set_pre_MJD0_ix(self.mjd0)

    def load(self, input_dir: str, control_index: int = 0):
        self.lcs[control_index] = AveragedLightCurve(
            self.colnames,
            control_index=control_index,
            filt=self.filt,
            mjdbinsize=self.mjdbinsize,
        )

        self.lcs[control_index].load_lc(input_dir, self.tnsname)

        if self.mjd0 is not None:
            self.set_pre_MJD0_ix(control_index=control_index)

    def load_all(self, input_dir: str, num_controls: int = 0):
        self.lcs = {}
        self.num_controls = 0

        print(
            f"\nLoading averaged SN light curve and {num_controls} averaged control light curves..."
        )

        # load averaged SN light curve
        self.load(input_dir)

        control_indices = find_all_control_indices(
            input_dir, self.tnsname, filt=self.filt
        )
        if len(control_indices) < num_controls:
            raise RuntimeError(
                f"Tried to load {num_controls} control light curves, but only {len(control_indices)} found: {control_indices}"
            )

        if num_controls > 0:
            # keep iterating over control indices until we successfully load num_controls averaged light curves
            control_index = 1
            while self.num_controls < num_controls:
                try:
                    self.load(input_dir, control_index=control_index)
                    self.num_controls += 1
                except:
                    print(
                        f"Could not load averaged control light curve {control_index}; skipping..."
                    )
                    del self.lcs[control_index]
                control_index += 1

        print(
            f"Successfully loaded averaged SN light curve and {self.num_controls} averaged control light curves (control indices: {self.control_lc_indices})"
        )

    def save_all(self, output_dir: str, overwrite: bool = False):
        print(
            f"\nDropping extra columns and saving averaged SN light curve and {self.num_controls} averaged control light curves..."
        )
        for control_index in self.lc_indices:
            self.lcs[control_index].drop_extra_columns()
            self.lcs[control_index].save_lc(
                output_dir, self.tnsname, overwrite=overwrite
            )
        print("Success")

    def __str__(self):
        return f"Averaged SN {self.tnsname} at {self.coords}: MJD0 = {self.mjd0}, {self.num_controls} control light curves"


# contains either o-band or c-band measurements only
class LightCurve(pdastrostatsclass):
    def __init__(
        self, colnames: PresetColumnNames, control_index=0, filt="o", **kwargs
    ):
        pdastrostatsclass.__init__(self, **kwargs)
        self.control_index = control_index
        self.filt = filt

        self.colnames = colnames
        self.colnames.add("dflux_new", self.colnames.dflux, overwrite=True)

    def set_df(self, t: pd.DataFrame):
        self.t = deepcopy(t)

    def get_preMJD0_indices(self, mjd0: float):
        return self.ix_inrange(
            colnames=self.colnames.mjd, uplim=mjd0, exclude_uplim=True
        )

    def get_postMJD0_indices(self, mjd0: float):
        return self.ix_inrange(colnames=self.colnames.mjd, lowlim=mjd0)

    def get_flags(self):
        return np.bitwise_or.reduce(self.t[self.colnames.mask])

    def get_good_indices(self, flag: Optional[int] = None):
        # if flag is 0, return all indices
        if flag == 0:
            return self.getindices()

        if flag is None:  # if no flag is given
            # check if mask column exists
            if not self.colnames.mask in self.t.columns:
                return self.getindices()

            # return all unmasked indices
            flag = self.get_flags()
        return self.ix_unmasked(self.colnames.mask, maskval=flag)

    def get_bad_indices(self, flag: Optional[int] = None):
        # if flag is 0, return no indices
        if flag == 0:
            return []

        if flag is None:  # if no flag is given
            # check if mask column exists
            if not self.colnames.mask in self.t.columns:
                return self.getindices()

            # return all masked indices
            flag = self.get_flags()
        return self.ix_masked(self.colnames.mask, maskval=flag)

    def get_ylims(
        self,
        indices: Optional[List[int]] = None,
        flag: Optional[int] = None,
        use_all: bool = False,
        mjd0: Optional[float] = None,
    ):
        if self.t.empty or len(self.t) < 2:
            return [None, None]

        if indices is None or len(indices) < 2:
            if use_all:
                indices = self.getindices()
            else:
                indices = self.get_good_indices(flag)
        if mjd0 is not None:
            pre_mjd0_ix = self.get_preMJD0_indices(mjd0)
            if len(pre_mjd0_ix) > 0:
                indices = AandB(indices, pre_mjd0_ix)

        flux_min = self.t.loc[indices, self.colnames.flux].min()
        flux_max = self.t.loc[indices, self.colnames.flux].max()
        offset = 0.05 * abs(flux_max - flux_min)

        return [flux_min - offset, flux_max + offset]

    def get_xlims(self, mjd0: Optional[float] = None, colname_attr: str = "mjd"):
        if self.t.empty or len(self.t) < 2:
            return [None, None]

        if colname_attr and not hasattr(self.colnames, colname_attr):
            raise AttributeError(f"Invalid colname_attr: '{colname_attr}'")
        colname = getattr(self.colnames, colname_attr)

        first_mjd = self.t.at[0, colname]
        last_mjd = self.t.at[len(self.t) - 1, colname]
        if mjd0 is None or mjd0 >= last_mjd or mjd0 <= first_mjd:
            end = last_mjd
        else:
            if colname_attr != "mjd" and colname_attr != "mjdbin":
                raise ValueError(
                    f"PresetColumnName attribute must be 'mjd' or 'mjdbin' if MJD0 is passed (received {colname_attr})"
                )
            end = mjd0

        return [first_mjd, end]

    def can_plot(self, ix: List[int], columns: List[str] = None):
        if columns is None:
            columns = [self.colnames.mjd, self.colnames.flux, self.colnames.dflux_new]

        available_columns = [col for col in columns if col in self.t.columns]
        if not available_columns:
            return False

        # check that we are plotting at least one row
        # and that the columns to plot are not all NaN values
        return len(ix) > 0 and not self.t.loc[ix, available_columns].isna().all().all()

    def remove_invalid_rows(self, verbose=False):
        dflux_zero_ix = self.ix_equal(colnames=[self.colnames.dflux], val=0)
        flux_nan_ix = self.ix_is_null(colnames=[self.colnames.flux])
        if len(AorB(dflux_zero_ix, flux_nan_ix)) > 0:
            if verbose:
                print(
                    f"Deleting {len(dflux_zero_ix) + len(flux_nan_ix)} rows with duJy=0 or uJy=NaN..."
                )
            self.t.drop(AorB(dflux_zero_ix, flux_nan_ix), inplace=True)

    def calculate_fdf_column(self, verbose=False):
        # replace infs with NaNs
        if verbose:
            print("Replacing infs with NaNs...")
        self.t.replace([np.inf, -np.inf], np.nan, inplace=True)

        # calculate flux/dflux
        if verbose:
            print(f"Calculating flux/dflux in for '{self.colnames.fdf}' column...")
        self.t[self.colnames.fdf] = (
            self.t[self.colnames.flux] / self.t[self.colnames.dflux_new]
        )

    def get_median_dflux(self, indices=None):
        if indices is None:
            indices = self.getindices()
        return np.nanmedian(self.t.loc[indices, self.colnames.dflux])

    def get_mean(
        self, colname: str, indices: List[int] = None, round_result: bool = False
    ) -> float:
        if indices is None:
            indices = self.getindices()

        self.calcaverage_sigmacutloop(
            colname, indices=indices, Nsigma=3.0, median_firstiteration=True
        )
        res = self.statparams["mean"]

        if res is None:
            print("WARNING: Could not converge on mean; taking median instead...")
            res = np.median(self.t.loc[indices, colname])

        if round_result:
            return round(res, 2)
        return res

    def get_stdev_flux(self, indices=None):
        self.calcaverage_sigmacutloop(
            self.colnames.flux, indices=indices, Nsigma=3.0, median_firstiteration=True
        )
        return self.statparams["stdev"]

    def add_noise_to_dflux(self, sigma_extra):
        new_dflux_colname = f"{self.colnames.dflux}_new"
        self.t[new_dflux_colname] = np.sqrt(
            self.t[self.colnames.dflux] * self.t[self.colnames.dflux] + sigma_extra**2
        )
        self.colnames.update("dflux_new", new_dflux_colname)
        self.calculate_fdf_column()

    def flag_by_control_stats(self, cut: ControlLightCurveCut):
        # flag SN measurements according to given bounds
        flag_x2_ix = self.ix_inrange(
            colnames=["c2_X2norm"], lowlim=cut.x2_max, exclude_lowlim=True
        )
        flag_stn_ix = self.ix_inrange(
            colnames=["c2_abs_stn"], lowlim=cut.snr_max, exclude_lowlim=True
        )
        flag_nclip_ix = self.ix_inrange(
            colnames=["c2_Nclip"], lowlim=cut.Nclip_max, exclude_lowlim=True
        )
        flag_ngood_ix = self.ix_inrange(
            colnames=["c2_Ngood"], uplim=cut.Ngood_min, exclude_uplim=True
        )
        self.update_mask_column(cut.x2_flag, flag_x2_ix)
        self.update_mask_column(cut.snr_flag, flag_stn_ix)
        self.update_mask_column(cut.Nclip_flag, flag_nclip_ix)
        self.update_mask_column(cut.Ngood_flag, flag_ngood_ix)

        # update mask column with control light curve cut on any measurements flagged according to given bounds
        zero_Nclip_ix = self.ix_equal("c2_Nclip", 0)
        unmasked_ix = self.ix_unmasked(
            self.colnames.mask,
            maskval=cut.x2_flag | cut.snr_flag | cut.Nclip_flag | cut.Ngood_flag,
        )
        self.update_mask_column(
            cut.questionable_flag, AnotB(unmasked_ix, zero_Nclip_ix)
        )
        self.update_mask_column(cut.flag, AnotB(self.getindices(), unmasked_ix))

    def copy_flags(self, flags_to_copy):
        self.t[self.colnames.mask] = self.t[self.colnames.mask].astype(np.int32)
        if len(self.t) < 1:
            return
        elif len(self.t) == 1:
            self.t.loc[0, self.colnames.mask] = (
                int(self.t.loc[0, self.colnames.mask]) | flags_to_copy
            )
        else:
            self.t[self.colnames.mask] = np.bitwise_or(
                self.t[self.colnames.mask], flags_to_copy
            )

    def average(
        self, cut: BadDayCut, previous_flags, mjdbinsize=1.0, flux2mag_sigmalimit=3.0
    ):
        avg_lc = AveragedLightCurve(
            self.colnames,
            self.control_index,
            filt=self.filt,
            mjdbinsize=mjdbinsize,
            columns=[
                self.colnames.mjd,
                self.colnames.mjdbin,
                self.colnames.flux,
                self.colnames.dflux,
                "stdev",
                "x2",
                "Nclip",
                "Ngood",
                "Nexcluded",
                self.colnames.mask,
            ],
            hexcols=[self.colnames.mask],
        )
        if self.control_index == 0:
            print(f"Now averaging SN light curve...")
        else:
            print(f"Now averaging control light curve {self.control_index}...")

        mjd = int(np.amin(self.t[self.colnames.mjd]))
        mjd_max = int(np.amax(self.t[self.colnames.mjd])) + 1

        while mjd <= mjd_max:
            range_ix = self.ix_inrange(
                colnames=[self.colnames.mjd],
                lowlim=mjd,
                uplim=mjd + mjdbinsize,
                exclude_uplim=True,
            )
            range_good_ix = self.ix_unmasked(
                self.colnames.mask, maskval=previous_flags, indices=range_ix
            )

            # add new row to averaged light curve
            new_row = {
                self.colnames.mjdbin: mjd + 0.5 * mjdbinsize,
                "Nclip": 0,
                "Ngood": 0,
                "Nexcluded": len(range_ix) - len(range_good_ix),
                self.colnames.mask: 0,
            }
            avglc_index = avg_lc.newrow(new_row)

            # if no measurements present, flag or skip over day
            if len(range_ix) < 1:
                avg_lc.update_mask_column(cut.flag, [avglc_index], remove_old=False)
                mjd += mjdbinsize
                continue

            # if no good measurements, average values anyway and flag
            if len(range_good_ix) < 1:
                # average flux
                self.calcaverage_sigmacutloop(
                    self.colnames.flux,
                    noisecol=self.colnames.dflux_new,
                    indices=range_ix,
                    Nsigma=3.0,
                    median_firstiteration=True,
                )
                fluxstatparams = deepcopy(self.statparams)

                # get average mjd
                self.calcaverage_sigmacutloop(
                    self.colnames.mjd,
                    indices=range_ix,
                    Nsigma=0,
                    median_firstiteration=False,
                )
                avg_mjd = self.statparams["mean"]

                # add row and flag
                row = {
                    self.colnames.mjd: avg_mjd,
                    self.colnames.flux: (
                        fluxstatparams["mean"]
                        if not fluxstatparams["mean"] is None
                        else np.nan
                    ),
                    self.colnames.dflux: (
                        fluxstatparams["mean_err"]
                        if not fluxstatparams["mean_err"] is None
                        else np.nan
                    ),
                    "stdev": (
                        fluxstatparams["stdev"]
                        if not fluxstatparams["stdev"] is None
                        else np.nan
                    ),
                    "x2": (
                        fluxstatparams["X2norm"]
                        if not fluxstatparams["X2norm"] is None
                        else np.nan
                    ),
                    "Nclip": (
                        fluxstatparams["Nclip"]
                        if not fluxstatparams["Nclip"] is None
                        else np.nan
                    ),
                    "Ngood": (
                        fluxstatparams["Ngood"]
                        if not fluxstatparams["Ngood"] is None
                        else np.nan
                    ),
                    self.colnames.mask: 0,
                }
                avg_lc.add2row(avglc_index, row)
                self.update_mask_column(cut.flag, range_ix, remove_old=False)
                avg_lc.update_mask_column(cut.flag, [avglc_index], remove_old=False)

                mjd += mjdbinsize
                continue

            # average good measurements
            self.calcaverage_sigmacutloop(
                self.colnames.flux,
                noisecol=self.colnames.dflux_new,
                indices=range_good_ix,
                Nsigma=3.0,
                median_firstiteration=True,
            )
            fluxstatparams = deepcopy(self.statparams)

            if fluxstatparams["mean"] is None or len(fluxstatparams["ix_good"]) < 1:
                self.update_mask_column(cut.flag, range_ix, remove_old=False)
                avg_lc.update_mask_column(cut.flag, [avglc_index], remove_old=False)
                mjd += mjdbinsize
                continue

            # get average mjd
            # TODO: SHOULD NOISECOL HERE BE DUJY OR NONE?
            self.calcaverage_sigmacutloop(
                self.colnames.mjd,
                noisecol=self.colnames.dflux_new,
                indices=fluxstatparams["ix_good"],
                Nsigma=0,
                median_firstiteration=False,
            )
            avg_mjd = self.statparams["mean"]

            # add row to averaged light curve
            row = {
                self.colnames.mjd: avg_mjd,
                self.colnames.flux: fluxstatparams["mean"],
                self.colnames.dflux: fluxstatparams["mean_err"],
                "stdev": fluxstatparams["stdev"],
                "x2": fluxstatparams["X2norm"],
                "Nclip": fluxstatparams["Nclip"],
                "Ngood": fluxstatparams["Ngood"],
                self.colnames.mask: 0,
            }
            avg_lc.add2row(avglc_index, row)

            # flag clipped measurements in lc
            if len(fluxstatparams["ix_clip"]) > 0:
                self.update_mask_column(
                    cut.ixclip_flag,
                    fluxstatparams["ix_clip"],
                    remove_old=False,
                )

            # if small number within this bin, flag measurements
            if len(range_good_ix) < 3:
                self.update_mask_column(cut.smallnum_flag, range_ix, remove_old=False)
                avg_lc.update_mask_column(
                    cut.smallnum_flag, [avglc_index], remove_old=False
                )
            # else check sigmacut bounds and flag
            else:
                is_bad = False
                if fluxstatparams["Ngood"] < cut.Ngood_min:
                    is_bad = True
                if fluxstatparams["Nclip"] > cut.Nclip_max:
                    is_bad = True
                if (
                    not (fluxstatparams["X2norm"] is None)
                    and fluxstatparams["X2norm"] > cut.x2_max
                ):
                    is_bad = True
                if is_bad:
                    self.update_mask_column(cut.flag, range_ix, remove_old=False)
                    avg_lc.update_mask_column(cut.flag, [avglc_index], remove_old=False)

            mjd += mjdbinsize

        avg_lc.flux2mag(
            self.colnames.flux,
            self.colnames.dflux,
            self.colnames.mag,
            self.colnames.dmag,
            zpt=23.9,
            upperlim_Nsigma=flux2mag_sigmalimit,
        )

        # TODO: not sure if needed
        for col in ["Nclip", "Ngood", "Nexcluded", self.colnames.mask]:
            avg_lc.t[col] = avg_lc.t[col].astype(np.int32)

        return avg_lc

    def apply_cut(self, cut: Cut):
        if not cut.can_apply_directly():
            raise RuntimeError(f"Cannot directly apply the following cut: {cut}")
        if not cut.column in self.t.columns:
            raise RuntimeError(
                f"No column name '{cut.column}' exists in light curve; cannot apply cut"
            )

        all_ix = self.getindices()
        kept_ix = self.ix_inrange(
            colnames=[cut.column], lowlim=cut.min_value, uplim=cut.max_value
        )
        cut_ix = AnotB(all_ix, kept_ix)

        self.update_mask_column(cut.flag, cut_ix)

        percent_cut = 100 * len(cut_ix) / len(all_ix)
        return percent_cut

    def remove_flag(self, flag):
        self.t[self.colnames.mask] = np.bitwise_and(
            self.t[self.colnames.mask].astype(int), ~flag
        )

    def update_mask_column(self, flag, indices, remove_old=True):
        if remove_old:
            # remove any old flags of the same value
            self.remove_flag(flag)

        if len(indices) > 1:
            flag_arr = np.full(self.t.loc[indices, self.colnames.mask].shape, flag)
            self.t.loc[indices, self.colnames.mask] = np.bitwise_or(
                self.t.loc[indices, self.colnames.mask].astype(int), flag_arr
            )
        elif len(indices) == 1:
            self.t.loc[indices, self.colnames.mask] = (
                int(self.t.loc[indices[0], self.colnames.mask]) | flag
            )

    def _clear_flux_offset_column(self):
        self.colnames.add("fluxoffset", f"{self.colnames.flux}_offset", overwrite=True)

        if self.colnames.fluxoffset in self.t.columns:
            print("Subtracting previous offset from flux column...")
            self.t[self.colnames.flux] -= self.t[self.colnames.fluxoffset]

        print("Setting current flux offset to 0...")
        self.t[self.colnames.fluxoffset] = 0

    def _get_region_mean(self, region_ix: List[int], maskval: int = None) -> float:
        indices = region_ix
        if not maskval is None and self.colnames.mask in self.t.columns:
            indices = self.ix_unmasked(
                self.colnames.mask, maskval=maskval, indices=region_ix
            )
        return self.get_mean(self.colnames.flux, indices=indices)

    def _add_offset(self, offset: int, region_ix: List[int]):
        offset_array = np.full(len(region_ix), offset)
        self.t.loc[region_ix, self.colnames.flux] += offset_array

        if not self.colnames.fluxoffset in self.t.columns:
            self.t[self.colnames.fluxoffset] = 0
        self.t.loc[region_ix, self.colnames.fluxoffset] += offset_array

    def _calculate_offset(
        self,
        ix1: List[int],
        ix2: List[int],
        num_measurements: int = 40,
        maskval: int = None,
    ) -> int:
        """Calculate the mean difference between two sets of measurements."""
        mean1 = self._get_region_mean(ix1[-num_measurements:], maskval=maskval)
        mean2 = self._get_region_mean(ix2[:num_measurements], maskval=maskval)
        return round(mean2 - mean1)

    def _get_region_indices(self):
        """Return indices for three regions based on MJD time intervals."""
        return {
            "1": self.ix_inrange(self.colnames.mjd, uplim=TEMPLATE_CHANGE_1_MJD),
            "2": self.ix_inrange(
                self.colnames.mjd,
                lowlim=TEMPLATE_CHANGE_1_MJD,
                uplim=TEMPLATE_CHANGE_2_MJD,
            ),
            "3": self.ix_inrange(self.colnames.mjd, lowlim=TEMPLATE_CHANGE_2_MJD),
        }

    def _get_offsets(
        self,
        mjd0: float,
        region_ix_dict: Dict,
        maskval: int = None,
        num_measurements: int = 40,
    ) -> Dict:
        if mjd0 > 57600:
            global_ix = (region_ix_dict["global"])[:num_measurements]
        else:
            global_ix = (region_ix_dict["global"])[-num_measurements:]

        return {
            "1": self._calculate_offset(
                region_ix_dict["1"],
                region_ix_dict["2"],
                maskval=maskval,
                num_measurements=num_measurements,
            ),
            "2": None,
            "3": self._calculate_offset(
                region_ix_dict["2"],
                region_ix_dict["3"],
                maskval=maskval,
                num_measurements=num_measurements,
            ),
            "global": -1 * round(self._get_region_mean(global_ix, maskval=maskval)),
        }

    def _apply_offsets(
        self,
        region_ix_dict: Dict,
        offset_dict: Dict,
    ):
        output = []
        for i in region_ix_dict.keys():
            region_ix = region_ix_dict[i]
            offset = offset_dict[i]
            if len(region_ix) > 0 and offset is not None:
                self._add_offset(offset, region_ix)
                output.append(f"Corrective flux {offset:0.2f} uJy added to region {i}")
        return output

    def _manual_template_correction(
        self,
        region1_offset: int = None,
        region2_offset: int = None,
        region3_offset: int = None,
    ):
        self._clear_flux_offset_column()

        region_ix_dict = self._get_region_indices()
        offset_dict: Dict[str : Optional[int]] = {
            "1": region1_offset,
            "2": region2_offset,
            "3": region3_offset,
        }

        return self._apply_offsets(region_ix_dict, offset_dict)

    def _auto_template_correction(
        self, mjd0: float, maskval: int = None, num_measurements: int = 40
    ):
        self._clear_flux_offset_column()

        region_ix_dict = self._get_region_indices()
        region_ix_dict["global"] = self.getindices()
        offset_dict = self._get_offsets(
            mjd0, region_ix_dict, maskval=maskval, num_measurements=num_measurements
        )

        return self._apply_offsets(region_ix_dict, offset_dict)

    def apply_template_correction(
        self,
        mjd0: float,
        maskval: int = None,
        region1_offset: Optional[int] = None,
        region2_offset: Optional[int] = None,
        region3_offset: Optional[int] = None,
        num_measurements: int = 40,
    ):
        self.colnames.add("fluxoffset", f"{self.colnames.flux}_offset", overwrite=True)
        if not (
            region1_offset is None and region2_offset is None and region3_offset is None
        ):
            print("Proceeding with manual template correction...")
            return self._manual_template_correction(
                region1_offset=region1_offset,
                region2_offset=region2_offset,
                region3_offset=region3_offset,
            )
        else:
            print("Proceeding with automatic template correction...")
            return self._auto_template_correction(
                mjd0, maskval=maskval, num_measurements=num_measurements
            )

    def drop_extra_columns(self, verbose=False):
        dropcols = []
        for col in [
            "Noffsetlc",
            self.colnames.fdf,
            "__tmp_SN",
            "SNR",
            "SNRsum",
            "SNRsumnorm",
            "SNRsim",
            "SNRsimsum",
            "c2_mean",
            "c2_mean_err",
            "c2_stdev",
            "c2_stdev_err",
            "c2_X2norm",
            "c2_Ngood",
            "c2_Nclip",
            "c2_Nmask",
            "c2_Nnan",
            "c2_abs_stn",
        ]:
            if col in self.t.columns:
                dropcols.append(col)
        for col in self.t.columns:
            if re.search("^c\d_", col):
                dropcols.append(col)

        if len(dropcols) > 0:
            if verbose:
                print(
                    f'Dropping extra columns ({f"control light curve {str(self.control_index)}" if self.control_index > 0 else "SN light curve"}): ',
                    dropcols,
                )
            self.t.drop(columns=dropcols, inplace=True)

    def check_column_names(self, required_column_names):
        if self.t is None:
            return

        for column_name in required_column_names:
            if not column_name in self.t.columns:
                raise RuntimeError(f"Missing required column: {column_name}")

    def load_lc(self, input_dir, tnsname, cleaned=False):
        filename = get_filename(
            input_dir, tnsname, self.filt, self.control_index, cleaned=cleaned
        )
        self.load_lc_by_filename(filename)

    def load_lc_by_filename(self, filename):
        self.load_spacesep(
            filename, delim_whitespace=True, hexcols=[self.colnames.mask]
        )
        self.check_column_names(
            required_column_names=self.colnames.get_required_column_names()
        )

    def save_lc(self, output_dir, tnsname, indices=None, overwrite=False, cleaned=True):
        filename = get_filename(
            output_dir, tnsname, self.filt, self.control_index, cleaned=cleaned
        )
        self.save_lc_by_filename(filename, indices=indices, overwrite=overwrite)

    def save_lc_by_filename(self, filename, indices=None, overwrite=False):
        self.write(
            filename=filename,
            indices=indices,
            overwrite=overwrite,
            hexcols=[self.colnames.mask],
        )

    def __str__(self):
        return self.t.to_string()


class AveragedLightCurve(LightCurve):
    def __init__(
        self,
        colnames: PresetColumnNames,
        control_index=0,
        filt="o",
        mjdbinsize=1.0,
        **kwargs,
    ):
        LightCurve.__init__(self, colnames, control_index, filt, **kwargs)
        self.mjdbinsize = mjdbinsize

        self._pre_mjd0_ix = None
        self._valid_mjd_ix = None

    @property
    def valid_mjd_ix(self):
        if self._valid_mjd_ix is None:
            raise RuntimeError(
                "Call self.set_valid_mjd_ix() first before accessing self.valid_mjd_ix"
            )
        if len(self._valid_mjd_ix) < 1:
            raise RuntimeError(
                f"No valid MJD indices found in light curve (control index {self.control_index})"
            )
        return self._valid_mjd_ix

    @property
    def pre_mjd0_ix(self):
        if self._pre_mjd0_ix is None:
            raise RuntimeError(
                "Call self.set_pre_MJD0_ix() first before accessing self.pre_mjd0_ix"
            )
        if len(self._pre_mjd0_ix) < 1:
            raise RuntimeError(
                f"No pre-MJD0 indices found in light curve (control index {self.control_index})"
            )
        return self._pre_mjd0_ix

    def has_pre_mjd0_ix(self):
        return self._pre_mjd0_ix is not None and len(self._pre_mjd0_ix) > 0

    def has_valid_mjd_ix(self):
        return self._valid_mjd_ix is not None and len(self._valid_mjd_ix) > 0

    def set_valid_mjd_ix(self, mjd_ranges: List[List[float]]) -> List[int]:
        def in_range(value, mjd_ranges):
            return any(r[0] <= value <= r[1] for r in mjd_ranges)

        self._valid_mjd_ix = self.t.index[
            self.t[self.colnames.mjdbin].apply(lambda x: in_range(x, mjd_ranges))
        ].tolist()

        if len(self._valid_mjd_ix) < 1:
            raise RuntimeError(
                f"No valid MJD indices found in light curve (control index {self.control_index}) for ranges: {mjd_ranges}"
            )

    def set_pre_MJD0_ix(self, mjd0: float):
        self._pre_mjd0_ix = self.ix_inrange(
            colnames=self.colnames.mjdbin, uplim=mjd0, exclude_uplim=True
        )

        if len(self._pre_mjd0_ix) < 1:
            print(
                f"WARNING: No pre-MJD0 indices found in light curve (control index {self.control_index})"
            )

    def get_min_and_max_mjd(self):
        if self.t.empty:
            raise RuntimeError("Light curve empty; cannot return min or max MJD ")
        return (
            np.floor(self.t[self.colnames.mjdbin].iloc[0]),
            np.ceil(self.t[self.colnames.mjdbin].iloc[-1]),
        )

    def get_xlims(self, mjd0: float = None, colname_attr: str = "mjdbin"):
        return super().get_xlims(mjd0=mjd0, colname_attr=colname_attr)

    def load_lc_by_filename(self, filename):
        self.load_spacesep(filename, delim_whitespace=True, hexcols=["Mask"])
        self.check_column_names(
            required_column_names=self.colnames.get_required_column_names(
                is_averaged=True
            )
        )

    def load_lc(self, input_dir, tnsname):
        filename = get_filename(
            input_dir, tnsname, self.filt, self.control_index, self.mjdbinsize
        )
        self.load_lc_by_filename(filename)

    def save_lc(self, output_dir, tnsname, indices=None, overwrite=False):
        filename = get_filename(
            output_dir, tnsname, self.filt, self.control_index, self.mjdbinsize
        )
        self.save_lc_by_filename(filename, indices=indices, overwrite=overwrite)


class LimCutsTable:
    def __init__(self, lc: LightCurve, snr_bound, indices=None):
        self.t = None

        self.lc = lc
        if indices is None:
            indices = self.lc.getindices()
        self.indices = indices

        self.good_ix, self.bad_ix = self.get_goodbad_indices(snr_bound)

    def get_goodbad_indices(self, snr_bound):
        if not self.lc.colnames.fdf in self.lc.t.columns:
            self.lc.calculate_fdf_column()

        good_ix = self.lc.ix_inrange(
            colnames=[self.lc.colnames.fdf],
            lowlim=-snr_bound,
            uplim=snr_bound,
            indices=self.indices,
        )
        bad_ix = AnotB(self.indices, good_ix)
        return good_ix, bad_ix

    def get_keptcut_indices(self, x2_max):
        kept_ix = self.lc.ix_inrange(
            colnames=self.lc.colnames.chisquare, uplim=x2_max, indices=self.indices
        )
        cut_ix = AnotB(self.indices, kept_ix)
        return kept_ix, cut_ix

    def calculate_row(self, x2_max, kept_ix=None, cut_ix=None):
        if kept_ix is None or cut_ix is None:
            kept_ix, cut_ix = self.get_keptcut_indices(x2_max)
        data = {
            "PSF Chi-Square Cut": x2_max,
            "N": len(self.indices),
            "Ngood": len(self.good_ix),
            "Nbad": len(self.bad_ix),
            "Nkept": len(kept_ix),
            "Ncut": len(cut_ix),
            "Ngood,kept": len(AandB(self.good_ix, kept_ix)),
            "Ngood,cut": len(AandB(self.good_ix, cut_ix)),
            "Nbad,kept": len(AandB(self.bad_ix, kept_ix)),
            "Nbad,cut": len(AandB(self.bad_ix, cut_ix)),
            "Pgood,kept": 100 * len(AandB(self.good_ix, kept_ix)) / len(self.indices),
            "Pgood,cut": 100 * len(AandB(self.good_ix, cut_ix)) / len(self.indices),
            "Pbad,kept": 100 * len(AandB(self.bad_ix, kept_ix)) / len(self.indices),
            "Pbad,cut": 100 * len(AandB(self.bad_ix, cut_ix)) / len(self.indices),
            "Ngood,kept/Ngood": 100
            * len(AandB(self.good_ix, kept_ix))
            / len(self.good_ix),
            "Ploss": 100 * len(AandB(self.good_ix, cut_ix)) / len(self.good_ix),
            "Pcontamination": 100 * len(AandB(self.bad_ix, kept_ix)) / len(kept_ix),
        }
        return data

    def calculate_table(self, cut_start, cut_stop, cut_step):
        print(
            f"Calculating loss and contamination for chi-square cuts from {cut_start} to {cut_stop}..."
        )

        self.t = pd.DataFrame(
            columns=[
                "PSF Chi-Square Cut",
                "N",
                "Ngood",
                "Nbad",
                "Nkept",
                "Ncut",
                "Ngood,kept",
                "Ngood,cut",
                "Nbad,kept",
                "Nbad,cut",
                "Pgood,kept",
                "Pgood,cut",
                "Pbad,kept",
                "Pbad,cut",
                "Ngood,kept/Ngood",
                "Ploss",
                "Pcontamination",
            ]
        )

        # for different x2 cuts decreasing from 50
        for cut in range(cut_start, cut_stop + 1, cut_step):
            kept_ix, cut_ix = self.get_keptcut_indices(cut)
            percent_kept = 100 * len(kept_ix) / len(self.indices)
            if percent_kept < 10:
                # less than 10% of measurements kept, so no chi-square cuts beyond this point are valid
                continue
            row = self.calculate_row(cut, kept_ix=kept_ix, cut_ix=cut_ix)
            self.t = new_row(self.t, row)


# will contain measurements from both filters (o-band and c-band)
class FullLightCurve:
    def __init__(
        self, control_index=0, ra: str = None, dec: str = None, mjd0: float = None
    ):
        self.t = None
        self.mjd0 = mjd0
        self.coords = Coordinates(ra, dec)
        self.control_index = control_index
        self.filts = None

    def get_tns_data(self, tnsname, api_key, tns_id, bot_name):
        if self.coords.is_incomplete() or self.mjd0 is None or np.isnan(self.mjd0):
            print("Querying TNS for RA, Dec, and discovery date...")
            json_data = query_tns(tnsname, api_key, tns_id, bot_name)
            if json_data is None:
                print(f"Skipping...")
                return

            if self.coords.is_empty():
                self.coords = get_tns_coords_from_json(json_data)
                print(f"Setting coordinates to TNS coordinates: {self.coords}")

            if self.mjd0 is None or np.isnan(self.mjd0):
                self.mjd0 = get_tns_mjd0_from_json(json_data)
                print(
                    f"Setting MJD0 to TNS discovery date minus {DISC_DATE_BUFFER}: {self.mjd0}"
                )

    # download the full light curve from ATLAS
    def download(self, headers, lookbacktime=None, max_mjd=None):
        if lookbacktime:
            min_mjd = float(Time.now().mjd - lookbacktime)
        else:
            min_mjd = 50000.0

        if not max_mjd:
            max_mjd = float(Time.now().mjd)

        print(
            f"Downloading ATLAS light curve at {self.coords} from {min_mjd} MJD to {max_mjd} MJD..."
        )

        if min_mjd > max_mjd:
            raise RuntimeError(f"max MJD {max_mjd} cannot be than min MJD {min_mjd}.")

        while True:
            try:
                result = query_atlas(
                    headers,
                    self.coords.ra.angle.degree,
                    self.coords.dec.angle.degree,
                    min_mjd,
                    max_mjd,
                )
                break
            except Exception as e:
                print("Exception caught: " + str(e))
                print("Trying again in 20 seconds! Waiting...")
                time.sleep(20)
                continue
        self.t = result

    def get_filt_lens(self):
        total_len = len(self.t)
        filt_lens = {
            "o": len(np.where(self.t["F"] == "o")[0]),
            "c": len(np.where(self.t["F"] == "c")[0]),
        }
        return total_len, filt_lens

    # divide the light curve by filter and save into separate files
    def save(self, colnames: PresetColumnNames, input_dir, tnsname, overwrite=False):
        if self.t is None:
            raise RuntimeError(
                "Cannot save light curve that hasn't been downloaded yet."
            )

        lc = LightCurve(colnames, control_index=self.control_index)
        lc.set_df(self.t)

        # sort data by mjd
        lc.t = lc.t.sort_values(by=["MJD"], ignore_index=True)

        # remove rows with duJy=0 or uJy=NaN
        dflux_zero_ix = lc.ix_equal(colnames=["duJy"], val=0)
        flux_nan_ix = lc.ix_is_null(colnames=["uJy"])
        if len(AorB(dflux_zero_ix, flux_nan_ix)) > 0:
            print(
                f"Deleting {len(dflux_zero_ix) + len(flux_nan_ix)} rows with duJy=0 or uJy=NaN..."
            )
            lc.t = lc.t.drop(AorB(dflux_zero_ix, flux_nan_ix))

        for filt in ["o", "c"]:
            filename = get_filename(
                input_dir, tnsname, filt=filt, control_index=self.control_index
            )
            indices = lc.ix_equal(colnames=["F"], val=filt)
            print(
                f"Saving downloaded light curve with filter {filt} (length {len(indices)}) at {filename}..."
            )
            lc.save_lc_by_filename(filename, indices=indices, overwrite=overwrite)

    def __str__(self):
        return f"Full light curve at {self.coords}: control ID = {self.control_index}, MJD0 = {self.mjd0}"


"""
ADD SIMULATIONS AND APPLY ROLLING SUM TO AVERAGED LIGHT CURVE 
"""


class Simulation(ABC):
    def __init__(self, model_name=None, **kwargs):
        """
        Initialize the Simulation object.
        """
        self.model_name = model_name
        self.brightness = None

    @abstractmethod
    def get_sim_flux(self, mjds, brightness, **kwargs):
        """
        Compute the simulated flux for the given MJDs and peak apparent magnitude.

        :param mjds: List or array of MJDs into which we inject the Simulation.
        :param brightness: Desired brightness (e.g., peak apparent magnitude or flux) of the simulation.
        :param kwargs: Additional Simulation parameters (e.g., sigma_sim=1.0 and time_peak_mjd=56780.5 for Gaussian)

        :return: An array of simulated flux values corresponding to the input MJDs.
        """
        pass

    def __str__(self):
        return f'Simulation with model name "{self.model_name}": brightness = {format_float(self.brightness)}'


class SimDetecSupernova(AveragedSupernova):
    def __init__(
        self,
        colnames: PresetColumnNames,
        tnsname: str = None,
        mjdbinsize: float = 1.0,
        mjd0: float | None = None,
        filt: str = "o",
        flag: int = 0x800000,
        **kwargs,
    ):
        AveragedSupernova.__init__(
            self,
            colnames,
            tnsname=tnsname,
            mjd0=mjd0,
            mjdbinsize=mjdbinsize,
            filt=filt,
            flag=flag,
            **kwargs,
        )
        self.lcs: Dict[int, SimDetecLightCurve] = {}

    def get_all_fom(self, sigma_kern: float):
        self.apply_rolling_sums(
            sigma_kern,
            valid_ix=self.has_valid_mjd_ix(),
            pre_mjd0_ix=self.has_pre_mjd0_ix(),
        )

        fom_list = []
        for control_index in self.control_lc_indices:
            fom = self.lcs[control_index].t.loc[
                self.lcs[control_index].valid_mjd_ix, self.colnames.snrsumnorm
            ]
            if not fom.empty:
                fom_list.append(fom)

        if fom_list:
            return pd.concat(fom_list, ignore_index=True)
        return pd.Series(dtype=float)

    def get_all_fom_dict(self, sigma_kerns: List[float]):
        print(f"Getting all control FOM for MJD ranges {self._mjd_ranges}...")
        if self._mjd_ranges is None:
            raise RuntimeError(f"Valid MJD ranges cannot be None")

        res = {}
        for sigma_kern in sigma_kerns:
            res[sigma_kern] = self.get_all_fom(sigma_kern)
        return res

    def get_prelim_fom_limit_ranges(
        self,
        sigma_kerns: List[float],
    ):
        print(
            f"Calculating preliminary valid FOM limit ranges for sigma_kerns {sigma_kerns}..."
        )
        res = {sigma_kern: [0.0] for sigma_kern in sigma_kerns}

        if self._mjd_ranges is None:
            raise RuntimeError(
                "Set self._mjd_ranges before calling self.get_prelim_fom_limit_ranges()"
            )

        all_fom_dict = self.get_all_fom_dict(sigma_kerns)

        for sigma_kern in sigma_kerns:
            max_fom = round(max(all_fom_dict[sigma_kern]) + 0.01, 2)
            res[sigma_kern].append(max_fom)

        print(f"Valid FOM limit ranges: {res}")
        return all_fom_dict, res

    def get_n_falsepos(
        self,
        sigma_kern: float,
        fom_limit: float,
        control_index: int = 0,
        verbose: bool = False,
    ):
        return self.lcs[control_index].get_n_falsepos(
            sigma_kern,
            fom_limit,
            self.mjd0,
            flag=self.flag,
            verbose=verbose,
        )

    def apply_rolling_sums(
        self,
        sigma_kern: float,
        valid_ix: bool = False,
        pre_mjd0_ix: bool = False,
    ):
        if valid_ix and not self.has_valid_mjd_ix():
            raise RuntimeError(
                "Valid MJD indices missing; set valid_ix=False or call self.set_valid_mjd_ix()"
            )
        if pre_mjd0_ix and not self.has_pre_mjd0_ix():
            raise RuntimeError(
                "Pre-MJD0 indices missing; set pre_mjd0_ix=False or call self.set_pre_MJD0_ix()"
            )

        msg = f"Applying rolling sum of sigma_kern={format_float(sigma_kern)} to all light curves"
        out = []
        sn_indices = self.lcs[0].getindices()
        if valid_ix:
            out.append("using only MJDs in included MJD ranges for all light curves")
            sn_indices = AandB(sn_indices, self.lcs[0].valid_mjd_ix)
        if pre_mjd0_ix:
            out.append("using only pre-MJD0 MJDs for SN light curve")
            sn_indices = AandB(sn_indices, self.lcs[0].pre_mjd0_ix)
        if out:
            msg += " (" + "; ".join(out) + ")"
        print(msg + "...")

        # apply rolling sum to SN lc
        self.lcs[0].apply_rolling_sum(sigma_kern, flag=self.flag, indices=sn_indices)

        # apply rolling sum to control lcs, filtering by valid MJD ranges if needed
        for control_index in self.control_lc_indices:
            self.lcs[control_index].apply_rolling_sum(
                sigma_kern,
                flag=self.flag,
                indices=self.lcs[control_index].valid_mjd_ix if valid_ix else None,
            )

    def remove_rolling_sums(self):
        for control_index in self.lc_indices:
            self.lcs[control_index].remove_rolling_sum()

    def remove_simulations(self):
        for control_index in self.lc_indices:
            self.lcs[control_index].remove_simulations()

    def load(self, input_dir: str, control_index: int = 0):
        self.lcs[control_index] = SimDetecLightCurve(
            self.colnames,
            control_index=control_index,
            filt=self.filt,
            mjdbinsize=self.mjdbinsize,
        )

        self.lcs[control_index].load_lc(input_dir, self.tnsname)

        if self.mjd0 is not None:
            self.set_pre_MJD0_ix(control_index=control_index)


class SimDetecLightCurve(AveragedLightCurve):
    def __init__(
        self,
        colnames: PresetColumnNames,
        control_index=0,
        filt="o",
        mjdbinsize=1.0,
        **kwargs,
    ):
        colnames.add_many(
            {
                "snr": "SNR",
                "snrsum": "SNR_sum",
                "snrsumnorm": "SNR_sumnorm",
                "fluxsim": f"{colnames.flux}_sim",
                "snrsim": "SNR_sim",
                "snrsimsum": "SNR_simsum",
            }
        )

        AveragedLightCurve.__init__(
            self, colnames, control_index, filt, mjdbinsize, **kwargs
        )

        self.cur_sigma_kern = None

    def get_n_falsepos(
        self,
        sigma_kern: float,
        fom_limit: float,
        mjd0: float,
        flag=0x800000,
        verbose=False,
    ):
        if self._pre_mjd0_ix is None:
            self.set_pre_MJD0_ix(mjd0)

        # for control light curves, loop through all valid indices
        # for the SN light curve, only loop through valid indices before MJD0
        indices = self.valid_mjd_ix if self.has_valid_mjd_ix() else self.getindices()
        if self.control_index == 0 and self.has_pre_mjd0_ix():
            indices = AandB(indices, self.pre_mjd0_ix)
        self.apply_rolling_sum(sigma_kern, flag=flag, indices=indices)

        # find any triggers above the FOM limit
        count = 0
        mjds = []
        above_lim = False
        for k in indices:
            if self.t.at[k, self.colnames.snrsumnorm] > fom_limit:
                if not above_lim:
                    mjds.append(self.t.at[k, self.colnames.mjdbin])
                    count += 1
                above_lim = True
            else:
                above_lim = False

        if verbose and len(mjds) > 0:
            print(
                f"sigma_kern {sigma_kern}, FOM limit {fom_limit:0.2f}, control index {self.control_index}: {count} trigger(s) at MJDs {mjds}"
            )
        return count

    def remove_columns(self, colnames: List[str]):
        dropcols = []
        for col in colnames:
            if col in self.t.columns:
                dropcols.append(col)
        if len(dropcols) > 0:
            self.t.drop(columns=dropcols, inplace=True)

    # remove rolling sum columns
    def remove_rolling_sum(self):
        self.cur_sigma_kern = None
        self.remove_columns(
            [
                "__tmp_SN",
                self.colnames.snr,
                self.colnames.snrsum,
                self.colnames.snrsumnorm,
            ]
        )

    # remove simulation columns
    def remove_simulations(self):
        self.remove_columns(
            [
                "__tmp_SN",
                self.colnames.fluxsim,
                self.colnames.snrsim,
                self.colnames.snrsimsum,
            ]
        )

    def get_new_gaussian_sigma(self, sigma_kern: float):
        """
        new_gaussian_sigma = round(sigma_kern / self.mjdbinsize)
        for TESS/other lcs:
        - if ratio > 3, leave it
        - if ratio < 3, round to 1 decimal place
        - if ratio < .3, don't round at all, or round to 5 decimal places
        """
        ratio = sigma_kern / self.mjdbinsize
        if ratio > 3:
            return round(ratio)
        elif ratio < 0.3:
            return round(ratio, 5)
        else:  # 0.3 < ratio < 3
            return round(ratio, 1)

    # apply a rolling sum to the light curve and add SNR, SNRsum, and SNRsumnorm columns
    def apply_rolling_sum(
        self, sigma_kern: float, indices=None, flag=0x800000, verbose=False
    ):
        if sigma_kern < self.mjdbinsize:
            raise ValueError(
                f"Cannot apply rolling sum with sigma_kern ({sigma_kern} days) less than MJD bin size ({self.mjdbinsize} days)"
            )

        if indices is None:
            indices = self.getindices()
        if len(indices) < 1:
            raise RuntimeError("not enough measurements to apply simulated gaussian")
        good_ix = AandB(indices, self.ix_unmasked(self.colnames.mask, flag))

        self.remove_rolling_sum()
        self.cur_sigma_kern = sigma_kern
        self.t.loc[indices, self.colnames.snr] = 0.0
        self.t.loc[good_ix, self.colnames.snr] = (
            self.t.loc[good_ix, self.colnames.flux]
            / self.t.loc[good_ix, self.colnames.dflux]
        )

        new_gaussian_sigma = self.get_new_gaussian_sigma(sigma_kern)
        windowsize = int(6 * new_gaussian_sigma)
        halfwindowsize = int(windowsize * 0.5) + 1
        if verbose:
            print(
                f"Sigma: {sigma_kern:0.2f} days; MJD bin size: {self.mjdbinsize:0.2f} days; sigma: {new_gaussian_sigma:0.2f} bins; window size: {windowsize} bins"
            )

        # calculate the rolling SNR sum
        l = len(self.t.loc[indices])
        dataindices = np.array(range(l) + np.full(l, halfwindowsize))
        temp = pd.Series(
            np.zeros(l + 2 * halfwindowsize), name=self.colnames.snr, dtype=np.float64
        )
        temp[dataindices] = self.t.loc[indices, self.colnames.snr]
        SNRsum = temp.rolling(windowsize, center=True, win_type="gaussian").sum(
            std=new_gaussian_sigma
        )
        self.t.loc[indices, self.colnames.snrsum] = list(SNRsum[dataindices])

        # normalize it
        norm_temp = pd.Series(
            np.zeros(l + 2 * halfwindowsize), name="norm", dtype=np.float64
        )
        norm_temp[np.array(range(l) + np.full(l, halfwindowsize))] = np.ones(l)
        norm_temp_sum = norm_temp.rolling(
            windowsize, center=True, win_type="gaussian"
        ).sum(std=new_gaussian_sigma)
        self.t.loc[indices, self.colnames.snrsumnorm] = list(
            SNRsum.loc[dataindices]
            / norm_temp_sum.loc[dataindices]
            * max(norm_temp_sum.loc[dataindices])
        )

    def add_sim_flux(
        self,
        good_ix: List[int],
        sim_flux,
        cur_sigma_kern: float = None,
        verbose: bool = False,
        remove_old: bool = True,
    ):
        """
        Add simulated flux to the light curve ("uJysim" column) and add "SNRsim" and "SNRsimsum" columns.

        :param good_ix: Unmasked/unflagged indices of the light curve.
        :param sim_flux: Array of simulated flux to add to the light curve
        :param cur_sigma_kern: The current kernel size of the rolling sum.
        :param remove_old: Remove any old simulations before adding the simulated flux.
        """
        if cur_sigma_kern is None:
            cur_sigma_kern = self.cur_sigma_kern
        if cur_sigma_kern is None:
            raise RuntimeError(
                "No current sigma kern passed as argument or stored during previously applied rolling sum."
            )

        if remove_old:
            self.remove_simulations()
            self.t.loc[good_ix, self.colnames.fluxsim] = self.t.loc[
                good_ix, self.colnames.flux
            ]
        self.t.loc[good_ix, self.colnames.fluxsim] += sim_flux

        # make sure all bad rows have SNRsim = 0.0 so they have no impact on the rolling SNRsum
        self.t[self.colnames.snrsim] = 0.0
        # include only simulated flux in the SNR
        self.t.loc[good_ix, self.colnames.snrsim] = (
            self.t.loc[good_ix, self.colnames.fluxsim]
            / self.t.loc[good_ix, self.colnames.dflux]
        )

        new_gaussian_sigma = round(cur_sigma_kern / self.mjdbinsize)
        windowsize = int(6 * new_gaussian_sigma)
        halfwindowsize = int(windowsize * 0.5) + 1
        if verbose:
            print(
                f"Sigma: {cur_sigma_kern:0.2f} days; MJD bin size: {self.mjdbinsize:0.2f} days; new sigma: {new_gaussian_sigma:0.2f} bins; window size: {windowsize} bins"
            )

        # calculate the rolling SNR sum for SNR with simulated flux
        l = len(self.t)
        dataindices = np.array(range(l) + np.full(l, halfwindowsize))
        temp = pd.Series(
            np.zeros(l + 2 * halfwindowsize),
            name=self.colnames.snrsim,
            dtype=np.float64,
        )
        temp[dataindices] = self.t[self.colnames.snrsim]
        SNRsimsum = temp.rolling(windowsize, center=True, win_type="gaussian").sum(
            std=new_gaussian_sigma
        )
        self.t[self.colnames.snrsimsum] = list(SNRsimsum.loc[dataindices])

    # get max FOM (for simulated FOM, column=SNRsimsum; else column=SNRsumnorm)
    # of measurements within the given indices
    def get_max_fom(self, indices: List[int] = None):
        if indices is None:
            indices = self.getindices()

        if self.colnames.snrsimsum in self.t.columns:
            colname = self.colnames.snrsimsum
        elif self.colnames.snrsumnorm in self.t.columns:
            colname = self.colnames.snrsumnorm
        else:
            raise RuntimeError(f"No FOM column found (columns: {self.t.columns})")

        max_fom_idx = self.t.loc[indices, colname].idxmax()

        max_fom_mjd = self.t.loc[max_fom_idx, self.colnames.mjdbin]
        max_fom = self.t.loc[max_fom_idx, colname]
        return max_fom_mjd, max_fom
