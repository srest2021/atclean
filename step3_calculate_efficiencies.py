#!/usr/bin/env python

from abc import ABC
import argparse
from collections import defaultdict
from configparser import ConfigParser
from copy import deepcopy
import itertools
import re
import sys
from typing import Any, Dict, List, Optional, Self
import numpy as np
import pandas as pd
from scipy import interpolate
from lightcurve import SimDetecSupernova
from pdastro import pdastrostatsclass
from step1_generate_sim_tables import (
    BRIGHTNESS_PARAM_PREFIX,
    TIME_PARAM_PREFIX,
    ListParam,
    ParamType,
    Params,
    find_prefix_in_list,
    parse_colname_info,
    parse_config_params,
)
from step2_generate_detec_tables import (
    NON_PARAM_COLNAMES,
    SimDetecTables,
    get_brightness_values_from_dir,
    get_detec_tables_output_dir,
    get_matching_ix,
    mjd_range_type,
)
from utils import (
    AandB,
    PresetColumnNames,
    SnInfoTable,
    format_float,
    get_allowed_presets,
    hexstring_to_int,
    load_config,
    load_json_config,
    make_dir_if_not_exists,
    new_row,
    print_progress_bar,
    validate_fom_limits,
)


class ContaminationTable:
    def __init__(self):
        self.is_prelim = False

    def calculate_row(
        self,
        sn: SimDetecSupernova,
        sigma_kern: float,
        fom_limit: float,
    ) -> Dict:
        row = {
            "sigma_kern": sigma_kern,
            "fom_limit": fom_limit,
            "n_falsepos": 0,
            "n_pos_controls": 0,
            "pct_pos_controls": np.nan,
        }

        for control_index in sn.lc_indices:
            n_falsepos = sn.get_n_falsepos(
                sigma_kern, fom_limit, control_index=control_index
            )
            row[f"n_falsepos_{control_index:02d}"] = n_falsepos

            if control_index > 0 and n_falsepos > 0:
                row["n_pos_controls"] += 1
                row["n_falsepos"] += n_falsepos

        row["pct_pos_controls"] = round(
            100 * row["n_pos_controls"] / sn.num_controls, 2
        )
        return row

    def construct_prelim_t(
        self,
        sn: SimDetecSupernova,
        # mjd_ranges: List[List[float]],
        sigma_kerns: List[float],
        fom_limits: Dict[int, List[float]],
    ):
        print(
            "Calculating preliminary contamination table for valid MJD ranges and preliminary FOM limit ranges..."
        )
        print(f"Using preliminary FOM limit ranges: {fom_limits}")

        if sn._mjd_ranges is None:
            raise RuntimeError(
                "Call sn.set_mjd_ranges() before calling self.construct_prelim_t()"
            )

        try:
            self.t = pd.DataFrame()
            for sigma_kern in sigma_kerns:
                for fom_limit in fom_limits[sigma_kern]:
                    row = self.calculate_row(sn, sigma_kern, fom_limit)
                    self.t = new_row(self.t, row)

            # number of false positives should always be 0 for min fom limits
            invalid_rows = self.t.iloc[1::2][self.t.iloc[1::2]["n_falsepos"] != 0]
            if not invalid_rows.empty:
                raise ValueError(
                    f"Invalid `n_falsepos` values found at row(s): {invalid_rows.index.tolist()}"
                )

            self.is_prelim = True
            print("Success")
        except Exception as e:
            self.t = None
            self.is_prelim = False
            raise RuntimeError(
                f"Could not construct preliminary contamination table: {str(e)}"
            )

    def get_initial_limits(
        self, index: int, prelim_fom_limit_ranges: Dict[int, List[float]]
    ):
        sigma_kern = self.t.loc[index, "sigma_kern"]
        lower_limit = round(self.t.loc[index, "fom_limit"], 2)
        upper_limit = round(self.t.loc[index + 1, "fom_limit"], 2)

        assert (
            lower_limit == prelim_fom_limit_ranges[sigma_kern][0]
            and upper_limit == prelim_fom_limit_ranges[sigma_kern][1]
        ), "Mismatch in preliminary FOM limits"

        return sigma_kern, lower_limit, upper_limit

    def calculate(
        self,
        sn: SimDetecSupernova,
        prelim_fom_limit_ranges: Dict[int, List[float]],
        sigma_kerns: List[float],
        target_value: int = 2,
        n_steps: int = 15,
        verbose: bool = False,
        convergence_threshold: float = 0.01,
    ) -> Dict[float, float]:
        """
        Calculate contamination metrics and refine FOM limits to achieve the target contamination level.

        :param sn: SimDetecSupernova object containing light curve data.
        :param prelim_fom_limit_ranges: Preliminary FOM limit ranges for each sigma_kern.
        :param sigma_kerns: List of rolling sum kernel sizes.
        :param tgt_value: Target number of positive control light curves.
        :param n_steps: Maximum number of iterations for refining FOM limits.
        :param verbose: Whether to print detailed logs.
        :param convergence_threshold: Threshold for convergence of FOM limits.
        :return: Dictionary of refined FOM limits for each sigma_kern.
        """

        if sn._mjd_ranges is None:
            raise RuntimeError(
                "Call sn.set_mjd_ranges() before calling self.calculate()"
            )

        if self.t is None or not self.is_prelim:
            self.construct_prelim_t(sn, sigma_kerns, prelim_fom_limit_ranges)
        if verbose:
            print("Preliminary contamination table: ")
            print(self.t.to_string())

        print(
            f"Refining FOM limits with up to {n_steps} iterations to achieve target contamination of {target_value} positive control light curves..."
        )
        fom_limits = {}
        i = 0
        while i < len(self.t) - 1:
            sigma_kern, lower_limit, upper_limit = self.get_initial_limits(
                i, prelim_fom_limit_ranges
            )
            if verbose:
                print(f"\n--- sigma_kern = {sigma_kern} ---")
                print(f"Preliminary FOM range: [{lower_limit}, {upper_limit}]")

            best_row = None
            for step in range(n_steps):
                new_fom_limit = round((upper_limit + lower_limit) / 2, 2)
                new_row = self.calculate_row(sn, sigma_kern, new_fom_limit)
                cur_value = new_row["n_pos_controls"]
                if verbose:
                    print(
                        f"Step {step + 1:>2}/{n_steps}: "
                        f"FOM limit={new_fom_limit:.2f}, "
                        f"Contamination={cur_value}, "
                        f"Bounds=[{lower_limit:.2f}, {upper_limit:.2f}]"
                    )

                if cur_value == target_value:
                    # this one satisfies the target contamination value — keep track of it
                    best_row = new_row

                # update limits based on the current contamination value
                if cur_value <= target_value:
                    # not enough many positives — tighten upper bound
                    upper_limit = new_fom_limit
                else:
                    # too many positives — tighten lower bound
                    lower_limit = new_fom_limit

                # check for convergence
                if round(abs(upper_limit - lower_limit), 2) <= convergence_threshold:
                    if verbose:
                        if best_row:
                            print(
                                f"→ Converged at FOM={best_row['fom_limit']} (exact match)"
                            )
                        else:
                            print(
                                f"→ Converged at FOM={new_fom_limit} (best guess), Contamination={cur_value}"
                            )
                    break

            # finalize the FOM limit for this sigma_kern
            if best_row is not None:
                final_row = best_row
            else:
                # fallback: use latest (if none were valid)
                final_row = new_row

            self.t.loc[i + 1, :] = final_row
            fom_limits[sigma_kern] = final_row["fom_limit"]
            if verbose:
                print(
                    f"✔ Final FOM limit for sigma_kern={sigma_kern}: "
                    f"{fom_limits[sigma_kern]:.2f} "
                    f"(Contamination={final_row['n_pos_controls']})"
                )
                print("-" * 22)

            i += 2

        # drop extra rows
        self.t = self.t.iloc[1::2].reset_index(drop=True)

        if verbose:
            print("\nFinal contamination table: ")
            print(self.__str__())

        return fom_limits

    def get_fom_limits_from_t(self):
        if self.t.empty:
            return {}

        if self.t["sigma_kern"].duplicated().any():
            fom_limits: Dict[int, List[float]] = defaultdict(list)
            for i in range(len(self.t)):
                fom_limits[self.t.loc[i, "sigma_kern"]].append(
                    self.t.loc[i, "fom_limit"]
                )
        else:
            fom_limits: Dict[int, float] = {}
            for i in range(len(self.t)):
                fom_limits[self.t.loc[i, "sigma_kern"]] = self.t.loc[i, "fom_limit"]
        return fom_limits

    def load(self, detec_tables_dir: str, prelim: bool = False):
        filename = f"{detec_tables_dir}/contamination{'_prelim' if prelim else ''}.txt"
        print(f"Loading contamination table at {filename}...")
        try:
            self.t = pd.read_table(filename, sep="\s+")
        except Exception as e:
            raise RuntimeError(
                f"Could not load efficiency table at {filename}: {str(e)}"
            )

    def save(self, detec_tables_dir: str, prelim: bool = False):
        make_dir_if_not_exists(detec_tables_dir)
        filename = f"{detec_tables_dir}/contamination{'_prelim' if prelim else ''}.txt"
        print(f"Saving contamination table as {filename}...")
        self.t.to_string(filename, index=False)

    def __str__(self):
        return self.t.to_string()


class EfficiencyTable(pdastrostatsclass):
    def __init__(
        self,
        sigma_kerns: List[float],
        params: Params,
        **kwargs,
    ):
        """
        Initialize an EfficiencyTable.

        :param sigma_kerns: List of detection algorithm kernel sizes.
        :param params: Collection of parameter names and possible values.
        """
        pdastrostatsclass.__init__(self, **kwargs)
        self.sigma_kerns: List[float] = sigma_kerns
        self._fom_limits: Optional[Dict[float, List[float]]] = None
        self.params: Params = params
        self.setup()

    def setup(self):
        """
        Set up the table columns for sigma_kerns, brightnesses, and other known parameter values.
        The time column name will be skipped when constructing columns.
        """
        # make sure brightness sorted for MagnitudeThresholdTable calculations later on
        self.params.brightness_param.values.sort()

        all_params = {
            "sigma_kern": self.sigma_kerns,
            self.params.brightness_param.name: self.params.brightness_param.values,
            **{param.name: param.values for param in self.params.other_params()},
        }

        print(f"\nSetting up efficiency table with columns: {all_params.keys()}")

        keys, values = zip(*(all_params).items())
        combinations = list(itertools.product(*values))
        self.t = pd.DataFrame(combinations, columns=keys)

        col_order = ["sigma_kern", self.params.brightness_param.name]
        col_order += [col for col in self.t.columns if col not in col_order]
        self.t = self.t[col_order]

    def set_fom_limits(
        self,
        fom_limits: (
            List[float]
            | List[List[float]]
            | Dict[float, float]
            | Dict[float, List[float]]
        ),
    ):
        self._fom_limits = validate_fom_limits(fom_limits, self.sigma_kerns)

    def get_params_at_index(self, index: int, skip_time_col: bool = False) -> Dict:
        """
        Get a dictionary of the parameter column-value pairs of the Simulation object at a certain row.
        Any known non-parameter column names (including brightness and, optionally, time) will be skipped.

        :param index: Index of the table from which to get the parameter column-value pairs.
        :param skip_time_col: Whether to exclude time columns (i.e., those starting with TIME_PARAM_PREFIX).
        """
        colnames = [
            col
            for col in self.t.columns
            if not (
                (skip_time_col and col.startswith(TIME_PARAM_PREFIX))
                or col.startswith(BRIGHTNESS_PARAM_PREFIX)
                or col in NON_PARAM_COLNAMES
                or re.search("^pct_detec_", col)
            )
        ]
        return dict(self.t.loc[index, colnames])

    def get_possible_values(self, column: str):
        """
        Return all unique values in the specified column.

        :param column: The name of the column for which to retrieve unique values.
        """
        if column not in self.t.columns:
            raise ValueError(f"Column '{column}' not found in the table.")
        return self.t[column].unique().tolist()

    def calculate_efficiencies(
        self,
        sd: SimDetecTables,
        fom_limits: (
            List[float]
            | List[List[float]]
            | Dict[float, float]
            | Dict[float, List[float]]
        ),
        progress_bar: bool = True,
        **kwargs,
    ):
        """
        For each row in the efficiency table, compute efficiencies for the FOM limits corresponding to the given sigma_kern.

        :param sd: SimDetecTables object that contains simulation information, max FOM, and other data needed to run the detection algorithm. Each SimDetecTable corresponds to one sigma_kern x peak_appmag combination.
        :param fom_limits: List or dictionary with sublists or single values as FOM limits.

        :param kwargs: Arbitrary number of pairs of column = value, column = range, or column = list of ranges.
        Example usage for columns A, B, C: self.get_efficiencies(sd, fom_limits, A=2, B=[5, 6], C=[[1, 2], [3, 4]])
        """
        self.set_fom_limits(fom_limits)

        l = len(self.t)
        print("Calculating efficiencies...")
        if progress_bar:
            print_progress_bar(0, l, prefix="Progress:", suffix="Complete", length=50)

        for i in range(l):
            sigma_kern = self.t.loc[i, "sigma_kern"]
            brightness = self.t.loc[i, self.params.brightness_param.name]

            if kwargs:
                params = kwargs
            else:
                params = self.get_params_at_index(i)

            for fom_limit in self._fom_limits[sigma_kern]:
                try:
                    efficiency = sd.get_efficiency(
                        sigma_kern, brightness, fom_limit, **params
                    )
                except Exception as e:
                    raise RuntimeError(
                        f"Could not calculate efficiency for sigma_kern={format_float(sigma_kern)}, brightness={format_float(brightness)}, fom_limit={format_float(fom_limit)}: {str(e)}"
                    )
                self.t.loc[i, f"pct_detec_{format_float(fom_limit)}"] = efficiency

            if progress_bar:
                print_progress_bar(
                    i + 1, l, prefix="Progress:", suffix="Complete", length=50
                )

        print("Success")
        print(self.__str__())

    def get_subset(
        self, fom_limits: Optional[List[float]] = None, **kwargs
    ) -> pd.DataFrame:
        """
        Get a subset of the table where the columns match the given values.

        :param fom_limits: List of FOM limits to get columns for. Set to None for all FOM limit columns.

        :param kwargs: Arbitrary number of pairs of column = value, column = range, or column = list of ranges.
        Example usage for columns A, B, C: self.get_subset([5.2, 7.8], sigma_kern=2, sigma_sim=[5, 6])

        :return: Efficiency of the rows that match the criteria.
        """
        colnames: List[str] = [
            "sigma_kern",
            self.params.brightness_param.name,
        ] + self.params.other_names()

        if not fom_limits is None:
            for col in self.t.columns:
                for fom_limit in fom_limits:
                    if re.search(f"^pct_detec_{format_float(fom_limit)}", col):
                        colnames.append(col)
        else:
            for col in self.t.columns:
                if re.search("^pct_detec_", col):
                    colnames.append(col)

        matching_ix = get_matching_ix(self, **kwargs)
        if len(matching_ix) > 0:
            return self.t.loc[matching_ix, colnames]
        else:
            return pd.DataFrame()

    def reset_table(self):
        """
        Remove any previously calculated efficiency columns.
        """
        for col in self.t.columns:
            if re.search("^pct_detec_", col):
                self.t.drop(col, axis=1, inplace=True)

    def _merge_fom_limits(
        self,
        other_fom_limits: Dict[float, List[float]],
    ):
        if self._fom_limits is None:
            self._fom_limits = deepcopy(other_fom_limits)
            return

        for other_key, other_values in other_fom_limits.items():
            if other_key in self._fom_limits:
                self._fom_limits[other_key].extend(other_values)
            else:
                self._fom_limits[other_key] = list(other_values)

        for key in self._fom_limits:
            self._fom_limits[key] = sorted(set(self._fom_limits[key]))

    def merge_tables(self, other: Self):
        """
        Add table content, sigma_kerns, and fom_limits from another EfficiencyTable.
        - Merges sigma_kerns with deduplication.
        - Validates and merges _fom_limits (deduplicated, sorted).
        - Concatenates DataFrame `t` if both are non-empty.
        """
        if not isinstance(other, EfficiencyTable):
            raise RuntimeError(
                f"Cannot merge EfficiencyTable with object type: {type(other)}"
            )

        self.sigma_kerns = list(sorted(set(self.sigma_kerns + other.sigma_kerns)))

        if other._fom_limits:
            self._merge_fom_limits(
                validate_fom_limits(other._fom_limits, other.sigma_kerns)
            )

        if other.t is not None and not other.t.empty:
            if self.t is not None and not self.t.empty:
                self.t = pd.concat([self.t, other.t], ignore_index=True)
            else:
                self.t = deepcopy(other.t)

    def load(self, detec_tables_dir: str, model_name: str):
        filename = f"{detec_tables_dir}/efficiencies_{model_name}.txt"
        print(f"Loading efficiency table at {filename}...")
        try:
            self.load_spacesep(filename, delim_whitespace=True)
        except Exception as e:
            raise RuntimeError(
                f"Could not load efficiency table at {filename}: {str(e)}"
            )

    def save(self, detec_tables_dir: str, model_name: str):
        make_dir_if_not_exists(detec_tables_dir)
        filename = f"{detec_tables_dir}/efficiencies_{model_name}.txt"
        print(f"Saving efficiency table as {filename}...")
        self.write(filename=filename, overwrite=True, index=False)

    def __str__(self):
        return self.t.to_string()


class MultipleRootsFound(Exception):
    def __init__(self, roots):
        self.roots = roots
        super().__init__(f"Multiple roots found: {roots}")


class MagnitudeThresholdTable:
    def __init__(self):
        """
        Initialize a MagnitudeThresholdTable.
        """
        self.all: pd.DataFrame = None
        self.best: pd.DataFrame = None

        self.sigma_kerns: List[float] = None
        self.select_param_name: str = None
        self.percents: List[float] = None

    def get_limits(self) -> tuple[float, float]:
        """
        Compute the lower and upper y-axis limits for plotting, based on all available
        magnitude threshold values across all specified percentiles.
        """
        if self.all is None or self.all.empty:
            raise RuntimeError("Table of all magnitude thresholds cannot be empty")
        if self.percents is None:
            raise RuntimeError("Percents cannot be None")

        y_cols = [f"mag_threshold_{format_float(p)}" for p in self.percents]
        y_vals = pd.concat([self.all[col].dropna() for col in y_cols])

        if not y_vals.empty:
            y_min = y_vals.min()
            y_max = y_vals.max()
            margin = 0.05 * (y_max - y_min) if y_max != y_min else 0.1
            ylim_lower = y_min - margin
            ylim_upper = y_max + margin
        else:
            ylim_lower = 0
            ylim_upper = 1

        return ylim_lower, ylim_upper

    def get_subset(
        self, sigma_kern: float, percent: float, drop_nans: bool = True
    ) -> tuple[pd.Series, pd.Series]:
        """
        Return x (select_param_name) and y (magnitude threshold) for a given sigma_kern and percent efficiency.
        Optionally drops rows with NaNs.
        """
        if self.all is None or self.all.empty or self.select_param_name is None:
            raise ValueError("self.all and self.select_param_name must be set")

        mag_threshold_colname = f"mag_threshold_{format_float(percent)}"
        if mag_threshold_colname not in self.all.columns:
            raise ValueError(
                f"Column '{mag_threshold_colname}' for {format_float(percent)}% efficiency does not exist"
            )

        subset = self.all[self.all["sigma_kern"] == sigma_kern][
            [self.select_param_name, mag_threshold_colname]
        ]
        if drop_nans:
            subset = subset.dropna()

        return (
            subset[self.select_param_name],
            subset[mag_threshold_colname],
        )

    def get_mag_threshold(self, x: pd.Series, y: pd.Series, percent: float):
        """
        Interpolates the detection function and solves for the x-value (brightness)
        that corresponds to the desired detection percent.

        :param x (pd.Series): Brightness parameter values.
        :param y (pd.Series): Detection percentages.
        :param percent (float): Target detection percent.

        Returns:
            float: Brightness value at which detection percentage reaches the target.
        """
        lx = x.to_list()
        ly = y.to_list()

        ly_reduced = np.array(ly) - percent
        freduced = interpolate.UnivariateSpline(lx, ly_reduced, s=0)
        roots = freduced.roots()
        if len(roots) > 1:
            raise MultipleRootsFound(roots)
            # return np.nan
        if len(roots) < 1:
            return np.nan
        return roots[0]

    def _generate_combinations(self, e: EfficiencyTable, select_param_name: str):
        for sigma_kern in e.sigma_kerns:
            for select_param_value in e.get_possible_values(select_param_name):
                for fom_limit in e._fom_limits[sigma_kern]:
                    yield sigma_kern, select_param_value, fom_limit

    def _calculate_row(
        self,
        e: EfficiencyTable,
        sigma_kern: float,
        select_param_name: str,
        select_param_value: Any,
        fom_limit: float,
        brightness_param_name: str,
        percents: List[float],
    ) -> Dict:
        subset = e.get_subset(
            sigma_kern=sigma_kern,
            **{select_param_name: select_param_value},
            fom_limits=[fom_limit],
        )
        row = {
            "sigma_kern": sigma_kern,
            select_param_name: select_param_value,
            "fom_limit": fom_limit,
        }
        for p in percents:
            colname = f"mag_threshold_{format_float(p)}"
            try:
                row[colname] = self.get_mag_threshold(
                    subset[brightness_param_name],
                    subset[f"pct_detec_{format_float(fom_limit)}"],
                    p,
                )
            except MultipleRootsFound as ex:
                print(
                    f"WARNING: Multiple roots found for sigma_kern={sigma_kern}, "
                    f"{select_param_name}={select_param_value}, fom_limit={fom_limit}, "
                    f"percent={p}. Roots: {ex.roots}"
                )
                # TODO: SHOULD THIS BE np.nan OR ex.roots[0]?
                row[colname] = ex.roots[0]
            except Exception as ex:
                print(
                    f"ERROR: Exception during mag threshold calc for sigma_kern={sigma_kern}, "
                    f"{select_param_name}={select_param_value}, fom_limit={fom_limit}, "
                    f"percent={p}. Exception: {ex}"
                )
                row[colname] = np.nan
        return row

    def _calculate_all(
        self,
        e: EfficiencyTable,
        select_param_name: str,
        percents: List[float] = [50, 80],
    ):
        """
        Compute the full magnitude threshold table across all combinations of sigma_kern,
        select_param, and FOM limit.
        """
        if e.t.empty:
            raise ValueError("EfficiencyTable cannot be empty")

        print("Calculating table of all magnitude thresholds...")
        columns = ["sigma_kern", select_param_name, "fom_limit"] + [
            f"mag_threshold_{format_float(p)}" for p in percents
        ]
        self.all = pd.DataFrame(columns=columns)

        brightness_param_name = find_prefix_in_list(
            e.t.columns, BRIGHTNESS_PARAM_PREFIX
        )

        for i, (sigma_kern, select_param_value, fom_limit) in enumerate(
            self._generate_combinations(e, select_param_name)
        ):
            row = self._calculate_row(
                e,
                sigma_kern,
                select_param_name,
                select_param_value,
                fom_limit,
                brightness_param_name,
                percents,
            )
            self.all.loc[i] = row
        print("Success")

    def _calculate_best(
        self,
        e: EfficiencyTable,
        select_param_name: str,
        percents: List[float] = [50, 80],
    ):
        """
        Extract the best (highest magnitude threshold) configuration for each select_param_value.
        """
        if e.t.empty:
            raise ValueError("EfficiencyTable cannot be empty")
        if self.all.empty:
            raise RuntimeError("Table of all magnitude thresholds cannot be empty")

        print("Calculating table of best magnitude thresholds...")
        columns = [select_param_name]
        for p in percents:
            columns.append(f"best_sigma_kern_{format_float(p)}")
            columns.append(f"best_mag_threshold_{format_float(p)}")
        self.best = pd.DataFrame(columns=columns)

        select_param_values = e.get_possible_values(select_param_name)

        for i in range(len(select_param_values)):
            select_param_value = select_param_values[i]
            select_param_value_ix = self.all[
                self.all[select_param_name] == select_param_value
            ].index

            best_data = {
                p: {"mag_threshold": -np.inf, "sigma_kern": np.nan} for p in percents
            }

            for sigma_kern in e.sigma_kerns:
                sigma_kern_ix = self.all[self.all["sigma_kern"] == sigma_kern].index
                ix = AandB(sigma_kern_ix, select_param_value_ix)

                for j in ix:
                    for p in percents:
                        m_key = f"mag_threshold_{format_float(p)}"
                        if self.all.loc[j, m_key] > best_data[p]["mag_threshold"]:
                            best_data[p]["mag_threshold"] = self.all.loc[j, m_key]
                            best_data[p]["sigma_kern"] = self.all.loc[j, "sigma_kern"]

            row = {select_param_name: select_param_value}
            for p in percents:
                row[f"best_mag_threshold_{format_float(p)}"] = best_data[p][
                    "mag_threshold"
                ]
                row[f"best_sigma_kern_{format_float(p)}"] = best_data[p]["sigma_kern"]
            self.best.loc[i] = row

        print("Success")

    def calculate(
        self, e: EfficiencyTable, select_param_name: str, percents: List[int] = [50, 80]
    ):
        """
        Orchestrates full and best-case magnitude threshold table calculations.

        :param e (EfficiencyTable): Source efficiency data.
        :param select_param_name (str): Name of the parameter to vary.
        :param percents (List[int]): Target detection percentages.
        """
        self.sigma_kerns = e.sigma_kerns
        self.select_param_name = select_param_name
        self.percents = percents

        self._calculate_all(e, select_param_name, percents)
        print(self.all.to_string(index=False))

        self._calculate_best(e, select_param_name, percents)
        print(self.best.to_string(index=False))

    def save(self, detec_tables_dir: str, model_name: str):
        make_dir_if_not_exists(detec_tables_dir)

        if self.all is not None and not self.all.empty:
            filename_all = (
                f"{detec_tables_dir}/all_magnitude_thresholds_{model_name}.txt"
            )
            print(f"Saving table of all magnitude thresholds as {filename_all}...")
            self.all.to_string(filename_all, index=False)

        if self.best is not None and not self.best.empty:
            filename_best = (
                f"{detec_tables_dir}/best_magnitude_thresholds_{model_name}.txt"
            )
            print(f"Saving table of best magnitude thresholds as {filename_best}...")
            self.best.to_string(filename_best, index=False)


class AnalysisLoop:
    def __init__(
        self, sigma_kerns: List[float], model_name: str, detec_tables_dir: str
    ):
        self.sigma_kerns = sigma_kerns
        self.model_name = model_name
        self.detec_tables_dir = detec_tables_dir

        self._sn: SimDetecSupernova = None
        self._tables: SimDetecTables = None
        self._params = Params()

        self.efficiencies: EfficiencyTable = None
        self.mag_thresholds: MagnitudeThresholdTable = None

    def _prepare_sn(
        self,
        mjd_ranges: Optional[List[List[float]]] = None,
        skip_control_ix: Optional[List] = None,
    ):
        if self._sn is None:
            raise RuntimeError(
                "Supernova (self._sn) must be set before calling self._prepare_sn()"
            )

        self._sn.remove_rolling_sums()
        self._sn.remove_simulations()
        if mjd_ranges is not None:
            self._sn.set_mjd_ranges(mjd_ranges)
        if skip_control_ix:
            print(f"Skipping control light curve indices: {skip_control_ix}")
            self._sn.remove_lc_indices(skip_control_ix)

    def load_sn(
        self,
        data_dir: str,
        colnames: PresetColumnNames,
        tnsname: str,
        num_controls: int,
        mjd0: Optional[float] = None,
        mjdbinsize: float = 1.0,
        filt: str = "o",
        mjd_ranges: Optional[List[List[float]]] = None,
        skip_control_ix: Optional[List] = None,
        flag: int = 0x800000,
    ):
        """
        Load the averaged SN and its control light curves.

        :param data_dir: Directory where the SN folder is located.
        :param tnsname: TNS name of the SN to load.
        :param num_controls: Number of averaged control light curves to load.
        :param mjd0: Discovery date or MJD of SN onset.
        :param mjdbinsize: MJD bin size of the averaged light curves to load.
        :param filt: Filter of the averaged light curves to load.
        :param mjd_ranges: Valid MJD ranges into which we inject Simulations.
        :param skip_control_ix: List of indices of control light curves which may NOT be randomly selected to have a Simulation injected.
        :param flag: Flag that denotes bad days in the binned light curves to load.
        """
        self._sn = SimDetecSupernova(
            colnames, tnsname, mjd0=mjd0, mjdbinsize=mjdbinsize, filt=filt, flag=flag
        )
        self._sn.load_all(data_dir, num_controls=num_controls)
        self._prepare_sn(mjd_ranges=mjd_ranges, skip_control_ix=skip_control_ix)

    def set_sn(
        self,
        sn: SimDetecSupernova,
        mjd_ranges: Optional[List[List[float]]] = None,
        skip_control_ix: Optional[List] = None,
    ):
        """
        Set a preloaded averaged SN and its control light curves.

        :param mjd_ranges: Valid MJD ranges into which we inject Simulations.
        :param skip_control_ix: List of indices of control light curves which may NOT be randomly selected to have a Simulation injected.
        """
        if not isinstance(sn, SimDetecSupernova):
            raise ValueError(
                "The provided object is not a valid SimDetecSupernova instance"
            )
        if sn.num_controls < 1 or len(sn.lcs) < 2:
            raise ValueError(
                "The SimDetecSupernova object must have at least one control light curve"
            )
        if not isinstance(sn.flag, int):
            raise ValueError(
                f"The flag attribute of the SimDetecSupernova object must be set to an integer (got {sn.flag})"
            )

        self._sn = deepcopy(sn)
        self._prepare_sn(mjd_ranges=mjd_ranges, skip_control_ix=skip_control_ix)

    def calculate_best_fom_limits(
        self, target_value: int = 2, n_steps: int = 15
    ) -> Dict[float, float]:
        if self._sn is None:
            raise RuntimeError(
                "Supernova (self._sn) must be set before calculating the best FOM limits"
            )
        if target_value < 1:
            raise ValueError("target_value must be >= 1")
        if n_steps < 2:
            raise ValueError("n_steps must be >= 2")

        print()
        _, prelim_fom_limit_ranges = self._sn.get_prelim_fom_limit_ranges(
            self.sigma_kerns
        )

        print()
        contam = ContaminationTable()
        contam.construct_prelim_t(self._sn, self.sigma_kerns, prelim_fom_limit_ranges)

        print()
        fom_limits = contam.calculate(
            self._sn,
            prelim_fom_limit_ranges,
            self.sigma_kerns,
            target_value=target_value,
            n_steps=n_steps,
            verbose=True,
        )
        contam.save(self.detec_tables_dir)

        print(f"Best FOM limits for each sigma_kern: {fom_limits}")
        return fom_limits

    def get_brightness_param_from_detec_tables(self, param_name: str = "brightness"):
        if self._sn is None:
            raise RuntimeError(
                "Supernova (self._sn) must be set before getting brightness parameter from SimDetecTables"
            )

        pattern = re.compile(
            rf"^simdetec_{re.escape(self.model_name)}_\d+\.\d+_(\d+\.\d+).({self._sn.filt})\.txt$"
        )
        values = get_brightness_values_from_dir(self.detec_tables_dir, pattern)
        if not values:
            raise FileNotFoundError(
                f"No matching detection tables found for brightness parameter (used model '{self.model_name}' and filter '{self._sn.filt}')."
            )

        brightness_param = ListParam(
            param_name, values, param_type=ParamType.BRIGHTNESS
        )
        self._params.add(brightness_param)
        print(self._params.brightness_param)

    def set_brightness_param(self, values: List[float], param_name="brightness"):
        if not values:
            raise ValueError("Brightness parameter values cannot be empty")

        brightness_param = ListParam(
            param_name, values, param_type=ParamType.BRIGHTNESS
        )
        self._params.add(brightness_param)

    def load_detec_tables_and_params(self):
        if self._sn is None:
            raise RuntimeError(
                "Supernova (self._sn) must be set before loading SimDetecTables"
            )
        if not self._params.has_brightness_param():
            raise RuntimeError(
                "Brightness parameter must be set before loading SimDetecTables"
            )

        self._tables = SimDetecTables(
            self._sn.filt,
            self._params.brightness_param,
            self.model_name,
            self.sigma_kerns,
        )
        self._tables.load_all(self.detec_tables_dir)
        self._params.merge(self._tables.get_params())

    def set_detec_tables(self, tables: SimDetecTables):
        self._tables = tables

    def set_params(self, params: Params):
        self._params = params
        self._params.validate(check_time=False)

    def calculate_efficiencies(
        self, target_value: int = 2, n_steps: int = 15
    ) -> EfficiencyTable:
        if self._params is None:
            raise RuntimeError(
                "Parameters (self._params) must be set before calculating efficiencies"
            )
        if self._tables is None:
            raise RuntimeError(
                "SimDetecTables (self._tables) must be set before calculating efficiencies"
            )

        fom_limits = self.calculate_best_fom_limits(
            target_value=target_value, n_steps=n_steps
        )
        self.efficiencies = EfficiencyTable(self.sigma_kerns, self._params)
        self.efficiencies.calculate_efficiencies(self._tables, fom_limits)
        self.efficiencies.save(self.detec_tables_dir, self.model_name)
        return self.efficiencies

    def calculate_mag_thresholds(
        self, select_param_name: str, percents: List[int] = [50, 80]
    ) -> MagnitudeThresholdTable:
        if self.efficiencies is None or self.efficiencies.t.empty:
            raise RuntimeError(
                "Efficiencies (self.efficiencies) must be calulated before calculating magnitude thresholds"
            )

        if self._params is None:
            raise RuntimeError(
                "Parameters (self._params) must be set before calculating magnitude thresholds"
            )
        if not self._params.has(select_param_name):
            raise ValueError(
                f"Parameter with name {select_param_name} is not known (known parameters: {self._params.all_names()})"
            )
        if (
            self._params.has_brightness_param()
            and self._params.brightness_param.name == select_param_name
        ):
            raise ValueError("Selected parameter name cannot be brightness")
        if (
            self._params.has_time_param()
            and self._params.time_param.name == select_param_name
        ):
            raise ValueError("Selected parameter name cannot be time")

        self.mag_thresholds = MagnitudeThresholdTable()
        self.mag_thresholds.calculate(
            self.efficiencies, select_param_name, percents=percents
        )
        self.mag_thresholds.save(self.detec_tables_dir, self.model_name)
        return self.mag_thresholds


class AtlasAnalysisLoop(AnalysisLoop):
    def __init__(self, sigma_kerns, model_name, detec_tables_dir):
        super().__init__(sigma_kerns, model_name, detec_tables_dir)

    def get_brightness_param_from_detec_tables(self):
        return super().get_brightness_param_from_detec_tables(param_name="peak_appmag")

    def calculate_mag_thresholds(
        self, select_param_name: Optional[str] = "sigma_sim", percents=[50, 80]
    ):
        return super().calculate_mag_thresholds(select_param_name, percents)


def define_args(
    config: ConfigParser, parser=None, usage=None, conflict_handler="resolve"
):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

    parser.add_argument("tnsname", type=str, help="transient name")
    parser.add_argument(
        "filter",
        type=str,
        default=None,
        help="filter of transient to inject simulations into",
    )
    parser.add_argument(
        "model_name", type=str, default="gaussian", help="name of model to use"
    )

    parser.add_argument(
        "--sigma_kerns",
        nargs="+",
        type=float,
        default=[5.0, 20.0, 40.0, 80.0, 100.0, 150.0, 200.0],
        help="list of kernel sizes in days for weighted gaussian rolling sum",
    )
    parser.add_argument(
        "--n_steps",
        type=int,
        default=15,
        help="maximum number of iterations for calculating best FOM limits via the bisection method",
    )
    parser.add_argument(
        "--n_pos_controls",
        type=int,
        default=2,
        help="target number of positive control light curves for a given FOM limit",
    )

    parser.add_argument(
        "-p",
        "--preset",
        type=str,
        default="atlas",
        help="preset name from config file (ex. atlas, rubin, tess)",
    )
    parser.add_argument(
        "--num_controls",
        type=int,
        default=int(config["download"]["num_controls"]),
        help="total number of averaged control light curves to load, not including skipped ones",
    )
    parser.add_argument(
        "--skip_control_ix",
        nargs="+",
        type=int,
        default=[],
        help="list of control indices to skip when loading control light curves",
    )
    parser.add_argument(
        "-m",
        "--mjd_bin_size",
        type=float,
        default=float(config["averaging"]["mjd_bin_size"]),
        help="MJD bin size in days of the target averaged light curves",
    )
    parser.add_argument(
        "--mjd0", type=float, default=None, help="transient start date in MJD"
    )
    parser.add_argument(
        "--mjd_ranges",
        type=mjd_range_type,
        default=None,
        help="List of MJD ranges as JSON string, e.g. '[[57233.5, 57328.5], [57466.5, 57535.5]]'",
    )

    return parser


if __name__ == "__main__":
    config = load_config("config.ini")
    args = define_args(config).parse_args()

    if " " in args.model_name:
        raise RuntimeError("Model name cannot have spaces.")

    allowed_presets = get_allowed_presets(config)
    if args.preset is None or args.preset not in allowed_presets:
        raise RuntimeError(
            f"Please specify the preset name to load from the config file (allowed presets: {allowed_presets})"
        )
    if args.filter in allowed_presets and args.filter != args.preset:
        print(
            f"WARNING: filter {args.filter} identified as preset in config.ini, but does not match arg preset {args.preset}"
        )
    print(f"\nLoading '{args.preset}' preset column names from config.ini...")
    colnames = PresetColumnNames(config, args.preset)
    print(colnames.__str__())

    mjd0 = args.mjd0
    if mjd0 is None:
        # get MJD0 from SnInfoTable
        sninfo = SnInfoTable(
            config["dir"]["output"], filename=config["dir"]["sninfo_filename"]
        )
        _, _, mjd0 = sninfo.get_info(args.tnsname)
    print(f"MJD0: {mjd0}")

    analysis_loop = AtlasAnalysisLoop(
        args.sigma_kerns,
        args.model_name,
        get_detec_tables_output_dir(config["dir"]["output"], args.tnsname),
    )

    analysis_loop.load_sn(
        config["dir"]["output"],
        colnames,
        args.tnsname,
        args.num_controls + len(args.skip_control_ix),
        mjd0=mjd0,
        mjdbinsize=float(args.mjd_bin_size),
        filt=args.filter,
        mjd_ranges=args.mjd_ranges,
        skip_control_ix=args.skip_control_ix,
        flag=hexstring_to_int(config["averaging"]["flag"]),
    )
    if args.mjd_ranges is not None:
        print(f"\nValid MJD ranges: {args.mjd_ranges}")

    print()
    analysis_loop.get_brightness_param_from_detec_tables()
    analysis_loop.load_detec_tables_and_params()
    analysis_loop.calculate_efficiencies(
        target_value=args.n_pos_controls, n_steps=args.n_steps
    )
    analysis_loop.calculate_mag_thresholds()
