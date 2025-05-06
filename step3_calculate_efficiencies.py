#!/usr/bin/env python

from collections import defaultdict
import itertools
import re
from typing import Dict, List, Optional, Self
import numpy as np
import pandas as pd
from scipy import interpolate
from lightcurve import SimDetecSupernova
from pdastro import pdastrostatsclass
from step1_generate_sim_tables import (
    BRIGHTNESS_PARAM_PREFIX,
    Params,
    find_prefix_in_list,
)
from step2_generate_detec_tables import SimDetecTables, get_matching_ix
from utils import (
    AandB,
    format_float,
    make_dir_if_not_exists,
    new_row,
    print_progress_bar,
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
        mjd_ranges: List[List[float]],
        sigma_kerns: List[float],
        fom_limits: Dict[int, List[float]],
    ):
        print(
            "Calculating preliminary contamination table for valid MJD ranges and preliminary FOM limit ranges..."
        )
        print(f"Using preliminary FOM limit ranges: {fom_limits}")

        try:
            sn.set_mjd_ranges(mjd_ranges)

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
        mjd_ranges: List[List[float]],
        sigma_kerns: List[float],
        target_value: int = 2,
        n_steps: int = 15,
        verbose: bool = False,
        convergence_threshold: float = 0.01,
    ):
        """
        Calculate contamination metrics and refine FOM limits to achieve the target contamination level.

        :param sn: SimDetecSupernova object containing light curve data.
        :param prelim_fom_limit_ranges: Preliminary FOM limit ranges for each sigma_kern.
        :param mjd_ranges: Valid MJD ranges for contamination calculation.
        :param sigma_kerns: List of rolling sum kernel sizes.
        :param tgt_value: Target number of positive control light curves.
        :param n_steps: Maximum number of iterations for refining FOM limits.
        :param verbose: Whether to print detailed logs.
        :param convergence_threshold: Threshold for convergence of FOM limits.
        :return: Dictionary of refined FOM limits for each sigma_kern.
        """

        if self.t is None or not self.is_prelim:
            self.construct_prelim_t(
                sn, mjd_ranges, sigma_kerns, prelim_fom_limit_ranges
            )
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
        fom_limits=None,
        **kwargs,
    ):
        """
        Initialize an EfficiencyTable.

        :param sigma_kerns: List of detection algorithm kernel sizes.
        :param params: Collection of parameter names and possible values.
        """
        pdastrostatsclass.__init__(self, **kwargs)
        self.sigma_kerns: List[float] = sigma_kerns
        self.fom_limits: Dict[float, List[float]] = self.set_fom_limits(fom_limits)
        self.params: Params = params

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

    def clear(self):
        self.t = None
        self.sigma_kerns = None
        self.fom_limits = None
        self.params = None

    # create dictionary of FOM limits, with sigma_kerns as the keys
    def validate_fom_limits(
        self,
        fom_limits: (
            List[float]
            | List[List[float]]
            | Dict[float, float]
            | Dict[float, List[float]]
        ),
    ) -> Dict[float, List[float]]:
        """
        Validate and convert FOM limits into a dictionary with sigma_kerns as keys.
        Supports:
        - List[float]: one FOM limit per sigma_kern
        - List[List[float]]: multiple FOM limits per sigma_kern
        - Dict[float, float]: one FOM limit per sigma_kern
        - Dict[float, List[float]]: multiple FOM limits per sigma_kern, already structured
        """
        # List[float] or List[List[float]]
        if isinstance(fom_limits, list):
            if not fom_limits:
                raise RuntimeError("No FOM limits provided")

            if len(fom_limits) != len(self.sigma_kerns):
                raise RuntimeError(
                    "Each entry in sigma_kerns must have a matching entry in fom_limits"
                )

            # List[float]
            if all(isinstance(x, (int, float)) for x in fom_limits):
                # wrap each float in a list
                return dict(zip(self.sigma_kerns, [[x] for x in fom_limits]))

            # List[List[float]]
            elif all(isinstance(x, list) for x in fom_limits):
                return dict(zip(self.sigma_kerns, fom_limits))

            else:
                raise TypeError("fom_limits list must contain sublists or numbers")

        # Dict[float, float] or Dict[float, List[float]]
        elif isinstance(fom_limits, dict):
            if set(fom_limits.keys()) != set(self.sigma_kerns):
                raise RuntimeError(
                    "FOM limits dict keys must exactly match sigma_kerns"
                )

            # wrap float values in lists if needed
            return {
                k: [v] if isinstance(v, (int, float)) else v
                for k, v in fom_limits.items()
            }

        else:
            raise TypeError(
                "fom_limits must be a list of lists/floats or a dict of lists/floats"
            )

    def set_fom_limits(
        self,
        fom_limits: (
            List[float]
            | List[List[float]]
            | Dict[float, float]
            | Dict[float, List[float]]
            | None
        ),
    ):
        if fom_limits is None:
            self.fom_limits = None
        else:
            self.fom_limits = self.validate_fom_limits(fom_limits)

    def get_params_at_index(self, index: int) -> Dict:
        """
        Get a dictionary of the parameter column-value pairs of the Simulation object at a certain row.
        Any known non-parameter column names (including brightness and time) will be skipped.

        :param index: Index of the table from which to get the parameter column-value pairs.
        """

        colnames = [col for col in self.t.columns if col in self.params.other_names()]
        return dict(self.t.loc[index, colnames])

    def get_possible_values(self, column: str):
        """
        Return all unique values in the specified column.

        :param column: The name of the column for which to retrieve unique values.
        """
        if column not in self.t.columns:
            raise ValueError(f"Column '{column}' not found in the table.")
        return self.t[column].unique().tolist()

    def get_efficiencies(
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

            for fom_limit in self.fom_limits[sigma_kern]:
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

    def merge_tables(self, other: Self):
        """
        Add table content, sigma_kerns, and fom_limits from another EfficiencyTable.
        """
        if not isinstance(other, EfficiencyTable):
            raise RuntimeError(
                f"Cannot merge EfficiencyTable with object type: {type(other)}"
            )

        self.sigma_kerns += other.sigma_kerns
        if not self.fom_limits is None:
            self.fom_limits.update(other.fom_limits)

        self.t = pd.concat([self.t, other.t], ignore_index=True)

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
            print(f"WARNING: Found more than one root {roots}; returning NaN")
            return np.nan
        if len(roots) < 1:
            return np.nan
        return roots[0]

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
        select_param_values = e.get_possible_values(select_param_name)

        i = 0
        for sigma_kern in e.sigma_kerns:
            for select_param_value in select_param_values:
                for fom_limit in e.fom_limits[sigma_kern]:
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
                        row[f"mag_threshold_{format_float(p)}"] = (
                            self.get_mag_threshold(
                                subset[brightness_param_name],
                                subset[f"pct_detec_{format_float(fom_limit)}"],
                                p,
                            )
                        )
                    self.all.loc[i] = row
                    i += 1

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
