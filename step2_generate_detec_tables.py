#!/usr/bin/env python

from abc import ABC, abstractmethod
from collections import defaultdict
from configparser import ConfigParser
import itertools
import os
import random
import argparse, re
from copy import deepcopy
import sys
from typing import Dict, List, Optional, Self, Tuple
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.optimize import root
from astropy.modeling.functional_models import Gaussian1D

from download import make_dir_if_not_exists
from pdastro import pdastrostatsclass
from step1_generate_sim_tables import (
    BRIGHTNESS_PARAM_PREFIX,
    TIME_PARAM_PREFIX,
    ListParam,
    Param,
    ParamType,
    Params,
    SimTable,
    SimTables,
    find_prefix_in_list,
    get_sim_tables_output_dir,
    parse_config_params,
    GAUSSIAN_MODEL_NAME,
    ASYMMETRIC_GAUSSIAN_MODEL_NAME,
)
from lightcurve import SimDetecLightCurve, SimDetecSupernova, Simulation
from utils import (
    AandB,
    PresetColumnNames,
    format_float,
    get_allowed_presets,
    hexstring_to_int,
    load_config,
    load_json_config,
    new_row,
    print_progress_bar,
)


NON_PARAM_COLNAMES = {
    "sigma_kern",
    "filename",
    "model_name",
    "mjd_colname",
    "mag_colname",
    "flux_colname",
    "control_index",
    "max_fom",
    "max_fom_mjd",
}


def get_detec_tables_output_dir(output_dir: str, tnsname: str):
    return os.path.join(output_dir, tnsname, "bump_analysis", "detec_tables")


# convert flux to magnitude
def flux2mag(flux: float):
    return -2.5 * np.log10(flux) + 23.9


# convert magnitude to flux
def mag2flux(mag: float):
    return 10 ** ((mag - 23.9) / -2.5)


def validate_kwarg_range(rng) -> bool:
    """
    Check if input is a valid [low, high] numeric range with low < high.
    """
    return (
        isinstance(rng, (list, tuple))
        and len(rng) == 2
        and all(isinstance(v, (int, float)) for v in rng)
        and rng[0] < rng[1]
    )


def get_matching_ix(table_object: pdastrostatsclass, **kwargs):
    """
    Get indices matching all conditions in kwargs.
    Supports:
    - Single value: A=3
    - Range: B=[1.0, 2.5] (must satisfy low < high)
    - List of ranges: C=[[1, 2], [3, 4]] (each must satisfy a < b)
    """
    matching_ix = table_object.getindices()
    for col, value in kwargs.items():
        if not col in table_object.t.columns:
            raise ValueError(
                f"Column '{col}' not found in table columns {table_object.t.columns}"
            )

        if isinstance(value, list):  # value is range or list of ranges
            # value is list of ranges
            if all(validate_kwarg_range(v) for v in value):
                mask = table_object.t.loc[matching_ix, col].apply(
                    lambda x: any(a <= x <= b for a, b in value)
                )
                matching_ix = table_object.t.index[mask].tolist()

            # value is range
            elif validate_kwarg_range(value):
                matching_ix = table_object.ix_inrange(
                    colnames=col,
                    lowlim=value[0],
                    uplim=value[1],
                    indices=matching_ix,
                )

            else:
                raise ValueError(f"Invalid range for column '{col}': {value}")

        else:  # single value
            matching_ix = table_object.ix_equal(col, value, indices=matching_ix)

    return matching_ix


class AsymmetricGaussian(Simulation):
    def __init__(self, model_name: str = ASYMMETRIC_GAUSSIAN_MODEL_NAME, **kwargs):
        """
        Initialize an asymmetric Gaussian simulation object.

        :param model_name: Name of the Gaussian model in the config file.
        """
        Simulation.__init__(self, model_name=model_name, **kwargs)
        self.g = None
        self.sigma_plus: float = None
        self.sigma_minus: float = None

    def new(self, sigma_plus: float, sigma_minus: float, peak_appmag: float):
        """
        :param sigma_plus: Sigma or kernel size of one half of the Gaussian.
        :param sigma_minus: Sigma or kernel size of the other half of the Gaussian.
        :param peak_appmag: Peak apparent magnitude of the Gaussian.
        """
        self.sigma_plus = sigma_plus
        self.sigma_minus = sigma_minus
        self.brightness = peak_appmag

        peak_flux = mag2flux(peak_appmag)
        x = np.arange(-100, 100, 0.01)
        g1 = Gaussian1D(amplitude=peak_flux, stddev=sigma_minus)(x)
        g2 = Gaussian1D(amplitude=peak_flux, stddev=sigma_plus)(x)

        ind = np.argmin(abs(x))
        g3 = np.copy(g1)
        g3[ind:] = g2[ind:]

        self.g = np.array([x, g3])

    def get_sim_flux(
        self,
        mjds,
        brightness: float,
        sigma_sim_plus: float = None,
        sigma_sim_minus: float = None,
        time_peak_mjd: float = None,
        **kwargs,
    ):
        """
        Get the interpolated function of the AsymmetricGaussian at a given peak MJD and match it to the given time array.

        :param mjds: Time array of MJDs.
        :param brightness: Desired peak apparent magnitude of the AsymmetricGaussian.
        :param sigma_plus: Sigma or kernel size of one half of the AsymmetricGaussian.
        :param sigma_minus: Sigma or kernel size of the other half of the AsymmetricGaussian.
        :param time_peak_mjd: MJD at which the AsymmetricGaussian should reach its peak apparent magnitude.

        :return: An array of simulated flux values corresponding to the input MJDs.
        """
        if sigma_sim_plus is None or sigma_sim_minus is None:
            raise RuntimeError(
                "sim_sigma_plus and sim_sigma_minus required to get flux of simulated asymmetric Gaussian."
            )
        if time_peak_mjd is None:
            raise RuntimeError(
                "Peak MJD required to get flux of simulated asymmetric Gaussian."
            )

        self.new(sigma_sim_plus, sigma_sim_minus, brightness)

        g = deepcopy(self.g)
        g[0, :] += time_peak_mjd

        fn = interp1d(g[0], g[1], bounds_error=False, fill_value=0)
        sim_flux = fn(mjds)
        return sim_flux

    def __str__(self):
        return (
            super().__str__()
            + f", sigma plus = {self.sigma_plus}, sigma_minus = {self.sigma_minus}"
        )


class Gaussian(AsymmetricGaussian):
    def __init__(self, model_name: str = GAUSSIAN_MODEL_NAME, **kwargs):
        """
        Initialize a Gaussian simulation object.

        :param model_name: Name of the Gaussian model in the config file.
        """
        AsymmetricGaussian.__init__(self, model_name=model_name, **kwargs)

    def get_sim_flux(
        self,
        mjds,
        brightness: float,
        sigma_sim: float = None,
        time_peak_mjd: float = None,
        **kwargs,
    ):
        """
        :param mjds: Time array of MJDs.
        :param brightness: Desired peak apparent magnitude of the Gaussian.
        :param sigma_sim: Desired sigma or kernel size of the Gaussian.
        :param time_peak_mjd: MJD at which the Gaussian should reach its peak apparent magnitude.
        """
        return super().get_sim_flux(
            mjds,
            brightness,
            sigma_sim_plus=sigma_sim,
            sigma_sim_minus=sigma_sim,
            time_peak_mjd=time_peak_mjd,
            **kwargs,
        )

    def __str__(self):
        return Simulation().__str__() + f", sigma = {self.sigma_plus}"


class Model(Simulation):
    def __init__(
        self,
        filename: str,
        mjd_colname: str = False,
        mag_colname: str = False,
        flux_colname: str = False,
        model_name: str = "pre_SN_outburst",
        **kwargs,
    ):
        """
        Initialize a model simulation object.

        :param filename: File name of the model to load.
        :param mjd_colname: MJD column name in the model file (None if present but no column name; False if not present).
        :param mag_colname: Magnitude column name in the model file (None if present but no column name; False if not present).
        :param flux_colname: Flux column name in the model file (None if present but no column name; False if not present).
        :param model_name: Name of the model assigned in the config file.
        """
        Simulation.__init__(self, model_name=model_name, **kwargs)
        self.t = None

        self.load(
            filename,
            mjd_colname=mjd_colname,
            mag_colname=mag_colname,
            flux_colname=flux_colname,
        )

    def load(
        self,
        filename: str,
        mjd_colname: str = False,
        mag_colname: str = False,
        flux_colname: str = False,
        verbose: bool = False,
    ):
        """
        Load the model from a file into a DataFrame.
        Discern which column is which using the given column names.
        If necessary, create any missing MJD, magnitude, or flux columns.

        :param filename: File name of the model to load.
        :param mjd_colname: MJD column name in the model file (None if present but no column name; False if not present).
        :param mag_colname: Magnitude column name in the model file (None if present but no column name; False if not present).
        :param flux_colname: Flux column name in the model file (None if present but no column name; False if not present).
        """
        if verbose:
            print(f"\nLoading model at {filename}...")

        if mag_colname is False and flux_colname is False:
            raise RuntimeError(
                f"Model must have either mag or flux column. Please set one or both fields to null or the correct column name."
            )

        try:
            header = "infer"
            if not mjd_colname and not mag_colname and not flux_colname:
                # all three column names are null or false
                header = None
            self.t = pd.read_table(filename, sep="\s+", header=header)
        except Exception as e:
            raise RuntimeError(f"Could not load model at {filename}: {str(e)}")

        if mjd_colname is False:
            # create MJD column and make it the first column
            columns = ["MJD"] + self.t.columns
            self.t["MJD"] = range(len(self.t))
            self.t = self.t[columns]
        else:
            # rename column to "MJD"
            self.t.rename(
                columns={0 if mjd_colname is None else mjd_colname: "MJD"}, inplace=True
            )

        if mag_colname is False:
            # flux column must be present
            # rename flux column to "uJy"
            self.t.rename(
                columns={1 if flux_colname is None else flux_colname: "uJy"},
                inplace=True,
            )
            # create mag column
            self.t["m"] = self.t["uJy"].apply(lambda flux: flux2mag(flux))
        else:
            # rename mag column to "m"
            self.t.rename(
                columns={1 if mag_colname is None else mag_colname: "m"}, inplace=True
            )
            if flux_colname is False:
                # create flux column
                self.t["uJy"] = self.t["m"].apply(lambda mag: mag2flux(mag))
            else:
                # rename flux column to "uJy"
                self.t.rename(
                    columns={2 if flux_colname is None else flux_colname: "uJy"},
                    inplace=True,
                )

        if verbose:
            print(self.t[["MJD", "m", "uJy"]].head().to_string())
            print("Success")

    def get_sim_flux(
        self, mjds, brightness: float, time_peak_mjd: float = None, **kwargs
    ):
        """
        Get the interpolated function of the model at a given peak MJD and peak apparent magnitude and match it to the given time array.

        :param mjds: Time array of MJDs.
        :param brightness: Desired peak apparent magnitude of the model.
        :param time_peak_mjd: MJD at which the model should reach its peak apparent magnitude.

        :return: The simulated flux array corresponding to the given time array.
        """
        self.brightness = brightness
        if time_peak_mjd is None:
            raise RuntimeError("Peak MJD required to construct simulated model.")

        # get original peak appmag index
        peak_idx = self.t["m"].idxmin()

        # scale flux to the desired peak appmag
        self.t["uJy"] *= mag2flux(brightness) / self.t.loc[peak_idx, "uJy"]

        # recalulate appmag column
        self.t["m"] = self.t["uJy"].apply(lambda flux: flux2mag(flux))

        # put peak appmag at peak_mjd
        self.t["MJD"] -= self.t.loc[peak_idx, "MJD"]
        self.t["MJD"] += time_peak_mjd

        # interpolate lc and match to time array
        fn = interp1d(self.t["MJD"], self.t["uJy"], bounds_error=False, fill_value=0)
        sim_flux = fn(mjds)
        return sim_flux

    def __str__(self):
        return super().__str__()


class SimDetecTable(SimTable):
    def __init__(self, sigma_kern: float, brightness: float, **kwargs):
        """
        Initialize a SimDetecTable.

        :sigma_kern: Sigma kernel of the algorithm.
        :brightness: Brightness (e.g., peak apparent magnitude or flux) for all simulations in this table.
        """
        SimTable.__init__(self, brightness, **kwargs)
        self.sigma_kern: float = sigma_kern

    def validate_model_name_col(self):
        if self.t.empty:
            print(
                "WARNING: Could not validate 'model_name' column because the SimDetecTable is empty"
            )
            return

        if not "model_name" in self.t.columns:
            print(
                "WARNING: Could not validate 'model_name' column because the SimDetecTable column does not exist"
            )
            return

        if self.t["model_name"].nunique() > 1:
            raise ValueError(
                "Multiple model types found in SimDetecTable, but expected only one"
            )

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
            )
        ]
        return dict(self.t.loc[index, colnames])

    def update_row_at_index(self, index: int, data: Dict):
        """
        Update a certain row of the table.

        :param index: Index of the table to update.
        :param data: Dictionary of column-value pairs.
        """
        for key, value in data.items():
            self.t.at[index, key] = value

    def get_detec_filename(self, model_name: str, detec_tables_dir: str) -> str:
        """
        Get the filename of the SimDetecTable.

        :param model_name: Name of the model of which the SimDetecTable contains simulations.
        :param detec_tables_dir: Directory where the SimDetecTable is located.
        """
        return f"{detec_tables_dir}/simdetec_{model_name}_{format_float(self.sigma_kern)}_{format_float(self.brightness)}.txt"

    def load_detec_table(self, model_name: str, detec_tables_dir: str):
        """
        Load an existing SimDetecTable.

        :param model_name: Name of the model of which the SimDetecTable contains simulations.
        :param detec_tables_dir: Directory where the SimDetecTable is located.
        """
        filename = self.get_detec_filename(model_name, detec_tables_dir)
        try:
            self.load_spacesep(filename, delim_whitespace=True)
        except Exception as e:
            raise RuntimeError(f"Could not load SimDetecTable at {filename}: {str(e)}")

    def load_from_sim_table(self, model_name: str, sim_tables_dir: str):
        """
        Load an existing SimTable and turn it into a SimDetecTable.

        :param model_name: Name of the model of which the SimTable contains simulations.
        :param sim_tables_dir: Directory where the SimTable is located.
        """
        super().load_sim_table(model_name, sim_tables_dir)
        self.t["sigma_kern"] = self.sigma_kern

    def save_detec_table(self, model_name: str, detec_tables_dir: str):
        """
        Save the current SimDetecTable.

        :param model_name: Name of the model of which the SimDetecTable contains simulations.
        :param detec_tables_dir: Directory where the SimDetecTable should be saved.
        """
        filename = self.get_detec_filename(model_name, detec_tables_dir)
        self.write(filename=filename, overwrite=True, index=False)

    def get_efficiency(self, fom_limit: float, **kwargs):
        """
        Get the efficiency where columns match all the given values and are within all the given ranges.

        :param kwargs: Arbitrary number of pairs of column = value, column = range, or column = list of ranges.
        Example usage for columns A, B, C: self.get_efficiency(10.0, A=2, B=[5, 6], C=[[1, 2], [3, 4]])

        :return: Efficiency of the rows that match the criteria.
        """
        if not "max_fom" in self.t.columns:
            raise ValueError("'max_fom' column not found in table")

        matching_ix = get_matching_ix(self, **kwargs)

        # no rows matched the params -> avoid division by 0
        if len(matching_ix) < 1:
            return 0.0

        detected_ix = self.ix_inrange("max_fom", lowlim=fom_limit, indices=matching_ix)
        efficiency = 100 * len(detected_ix) / len(matching_ix)
        return efficiency


class SimDetecTables:
    def __init__(
        self, brightness_param: Param, model_name: str, sigma_kerns: List[float]
    ):
        self.model_name: str = model_name
        self.sigma_kerns: List[float] = sigma_kerns
        self.brightness_param = brightness_param
        self.d: Dict[float, Dict[float, SimDetecTable]] = {}

    def get_table(self, sigma_kern: float, brightness: float):
        return self.d[sigma_kern][brightness]

    def update_row_at_index(
        self, sigma_kern: float, brightness: float, index: int, data: Dict
    ):
        self.d[sigma_kern][brightness].update_row_at_index(index, data)

    def get_efficiency(
        self, sigma_kern: float, brightness: float, fom_limit: float, **params
    ):
        return self.d[sigma_kern][brightness].get_efficiency(fom_limit, **params)

    def save_detec_table(
        self, sigma_kern: float, brightness: float, detec_tables_dir: str
    ):
        self.d[sigma_kern][brightness].save_detec_table(
            self.model_name, detec_tables_dir
        )

    def save_all(self, detec_tables_dir: str):
        print(f"\nSaving SimDetecTables in directory: {detec_tables_dir}")
        make_dir_if_not_exists(detec_tables_dir)
        for sigma_kern in self.d.keys():
            for table in self.d[sigma_kern].values():
                table.save_detec_table(self.model_name, detec_tables_dir)
        print("Success")

    def load_all_from_sim_tables(self, sim_tables_dir: str):
        """
        Load existing SimTables and turn them into SimDetecTables.

        :param sim_tables_dir: Directory where the SimTables are located.
        """
        print(
            f"\nConstructing SimDetecTables from existing SimTables in directory: {sim_tables_dir}"
        )

        if self.brightness_param is None:
            raise RuntimeError("brightness_param cannot be None")

        for sigma_kern in self.sigma_kerns:
            self.d[sigma_kern] = {}
            for brightness in self.brightness_param.values:
                self.d[sigma_kern][brightness] = SimDetecTable(sigma_kern, brightness)
                self.d[sigma_kern][brightness].load_from_sim_table(
                    self.model_name, sim_tables_dir
                )
        print("Success")

    def load_all(self, detec_tables_dir: str):
        """
        Load existing SimDetecTables.

        :param detec_tables_dir: Directory where the SimDetecTables are located.
        """
        print(f"\nLoading SimDetecTables from directory: {detec_tables_dir}")

        if self.brightness_param is None:
            raise RuntimeError("brightness_param cannot be None")

        for sigma_kern in self.sigma_kerns:
            self.d[sigma_kern] = {}
            for brightness in self.brightness_param.values:
                self.d[sigma_kern][brightness] = SimDetecTable(sigma_kern, brightness)
                self.d[sigma_kern][brightness].load_detec_table(
                    self.model_name, detec_tables_dir
                )
        print("Success")


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

    def get_limits(self):
        if self.all.empty:
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

    def get_mag_threshold(self, x: pd.Series, y: pd.Series, percent: float):
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


class SimDetecLoop(ABC):
    def __init__(self, sigma_kerns: List[float], **kwargs):
        self.sigma_kerns: List[float] = sigma_kerns
        self.brightness_param: Param = None

        self.sn: SimDetecSupernova = None
        self.e: EfficiencyTable = None
        self.sd: SimDetecTables = None

    def _get_brightness_values_from_dir(
        self,
        directory: str,
        pattern: re.Pattern,
    ):
        filenames = os.listdir(directory)

        brightnesses = set()

        for filename in filenames:
            match = pattern.match(filename)
            if match:
                brightness = float(match.group(1))
                brightnesses.add(brightness)

        res = list(brightnesses)
        res.sort()
        return res

    def get_brightness_param_from_detec_tables(
        self,
        model_name: str,
        detec_tables_dir: str,
        param_name: str = "brightness",
    ):
        pattern = re.compile(
            rf"^simdetec_{re.escape(model_name)}_\d+\.\d+_(\d+\.\d+)\.txt$"
        )
        values = self._get_brightness_values_from_dir(detec_tables_dir, pattern)
        self.brightness_param = ListParam(
            param_name, values, param_type=ParamType.BRIGHTNESS
        )
        print(self.brightness_param)

    def get_brightness_param_from_sim_tables(
        self,
        model_name: str,
        sim_tables_dir: str,
        param_name: str = "brightness",
    ):
        pattern = re.compile(rf"^sim_{re.escape(model_name)}_(\d+\.\d+)\.txt$")
        values = self._get_brightness_values_from_dir(sim_tables_dir, pattern)
        self.brightness_param = ListParam(
            param_name, values, param_type=ParamType.BRIGHTNESS
        )
        print(self.brightness_param)

    def set_brightness_param(self, values: List[float]):
        self.brightness_param = ListParam(
            "brightness", values, param_type=ParamType.BRIGHTNESS
        )

    def load_sn(
        self,
        data_dir: str,
        colnames: PresetColumnNames,
        tnsname: str,
        num_controls: int,
        mjdbinsize: float = 1.0,
        filt: str = "o",
        flag: int = 0x800000,
    ):
        """
        Load the averaged SN and its control light curves.

        :param data_dir: Directory where the SN folder is located.
        :param tnsname: TNS name of the SN to load.
        :param num_controls: Number of averaged control light curves to load.
        :param mjdbinsize: MJD bin size of the averaged light curves to load.
        :param filt: Filter of the averaged light curves to load.
        :param flag: Flag that denotes bad days in the binned light curves to load.
        """
        self.sn = SimDetecSupernova(
            colnames, tnsname, mjdbinsize=mjdbinsize, filt=filt, flag=flag
        )
        self.sn.load_all(data_dir, num_controls=num_controls)
        self.sn.remove_rolling_sums()
        self.sn.remove_simulations()

    def set_sn(self, sn: SimDetecSupernova):
        if not isinstance(sn, SimDetecSupernova):
            raise ValueError(
                "The provided object is not a valid SimDetecSupernova instance"
            )
        if sn.num_controls < 1 or len(sn.lcs) < 2:
            raise ValueError(
                "The SimDetecSupernova object must have at least one control light curve"
            )

        self.sn = sn
        self.sn.remove_rolling_sums()
        self.sn.remove_simulations()

    def load_sim_tables(self, model_name: str, sim_tables_dir: str):
        """
        Load existing SimTables and construct SimDetecTables out of them.

        :param model_name: Name of the model whose tables will be loaded.
        :param sim_tables_dir: Directory where the SimTables are located.
        """
        self.sd = SimDetecTables(self.brightness_param, model_name, self.sigma_kerns)
        self.sd.load_all_from_sim_tables(sim_tables_dir)

    def load_detec_tables(self, model_name: str, detec_tables_dir: str):
        """
        Load existing SimDetecTables.

        :param model_name: Name of the model whose tables will be loaded.
        :param detec_tables_dir: Directory where the SimDetecTables are located.
        """
        self.sd = SimDetecTables(self.brightness_param, model_name, self.sigma_kerns)
        self.sd.load_all(detec_tables_dir)

    def load_sim(self, table_row: Dict, verbose: bool = False) -> Simulation:
        """
        Construct and return a Simulation object given a row from a SimTable or SimDetecTable.
        """
        if verbose:
            print(
                "\tLoading Simulation from model info in first row of loaded SimTable"
            )
        model_name = table_row["model_name"]
        filename = (
            None
            if not isinstance(table_row["filename"], str)
            else table_row["filename"]
        )

        def get_col_val(colname, table_row):
            if colname in table_row:
                return None if np.isnan(table_row[colname]) else table_row[colname]
            else:
                return False

        mjd_colname = get_col_val("mjd_colname", table_row)
        mag_colname = get_col_val("mag_colname", table_row)
        flux_colname = get_col_val("flux_colname", table_row)

        if model_name == GAUSSIAN_MODEL_NAME:
            if verbose:
                print("\tUsing Gaussian simulations")
            sim = Gaussian()
        elif model_name == ASYMMETRIC_GAUSSIAN_MODEL_NAME:
            if verbose:
                print("\tUsing AsymmetricGaussian simulations")
            sim = AsymmetricGaussian()
        else:
            if verbose:
                print(
                    f"\tUsing '{model_name}' simulations with MJD column {mjd_colname}, mag column {mag_colname}, and flux column {flux_colname} at filename: {filename}"
                )
            sim = Model(
                filename=filename,
                mjd_colname=mjd_colname,
                mag_colname=mag_colname,
                flux_colname=flux_colname,
                model_name=model_name,
            )
        return sim

    def add_simulation_to_lc(
        self,
        sigma_kern: float,
        brightness: float,
        control_index: int,
        sim: Simulation,
        remove_old: bool = True,
        verbose: bool = False,
        **kwargs,
    ):
        """
        Add any Simulation object to a copy of a light curve, specifying parameters using keyword arguments.

        :param sigma_kern: The current sigma of the rolling sum.
        :param brightness: The desired brightness of the Simulation to add.
        :control_index: The control index of the light curve to add the Simulation to.
        :param sim: The Simulation to add.
        :param remove_old: Remove any old simulations before adding the new simulated flux.
        """
        if verbose:
            print(f"Adding simulation: {sim}")
            if kwargs:
                print(f"Additional simulation parameters: {kwargs}")

        lc = deepcopy(self.sn.lcs[control_index])
        good_ix = lc.get_good_indices(flag=self.sn.flag)

        sim_flux = sim.get_sim_flux(
            lc.t.loc[good_ix, self.sn.colnames.mjd], brightness, **kwargs
        )

        lc.add_sim_flux(
            good_ix,
            sim_flux,
            cur_sigma_kern=sigma_kern,
            verbose=verbose,
            remove_old=remove_old,
        )
        return lc

    @abstractmethod
    def get_max_fom_indices(
        self,
        sim_lc: SimDetecLightCurve,
        **kwargs,
    ):
        """
        From a light curve with a Simulation injected, return the indices within which to search for the max FOM.

        :param sim_lc: SimDetecLightCurve with a Simulation injected.
        """
        pass

    def update_sd_row(
        self,
        sigma_kern: float,
        brightness: float,
        index: int,
        control_index: int,
        max_fom: float,
        max_fom_mjd: float,
    ):
        """
        Update a certain row of a SimDetecTable with info about an injected Simulation, i.e., the index of the control light curve it was added to, the max FOM, and the max FOM MJD.

        :param sigma_kern: Sigma of the desired SimDetecTable.
        :param peak_appmag: Peak apparent magnitude of the desired SimDetecTable.
        :param index: Index of the row to update.
        :param control_index: Control light curve index to which the Simulation was added.
        :param max_fom: Max FOM of the simulated flux.
        :param max_fom_mjd: MJD of the max FOM of the simulated flux.
        """
        data = {
            "control_index": control_index,
            "max_fom": max_fom,
            "max_fom_mjd": max_fom_mjd,
        }
        self.sd.update_row_at_index(sigma_kern, brightness, index, data)

    def calculate_efficiencies(
        self,
        fom_limits: (
            List[float]
            | List[List[float]]
            | Dict[float, float]
            | Dict[float, List[float]]
        ),
        params: Params,
        detec_tables_dir: str,
        model_name: str,
        progress_bar: bool = True,
    ):
        """
        Construct and save an EfficiencyTable that contains efficiencies for every combination of a Simulation's sigma_kern, peak_appmag, and other parameters EXCEPT the time parameter.

        :param fom_limits: Dict or List of FOM limits.
        :param params: Collection of parameter names and possible values.
        :param detec_tables_dir: Directory where the EfficiencyTable should be saved.
        :param model_name: Name of the model for which to calculate efficiencies.
        """
        self.e = EfficiencyTable(self.sigma_kerns, params)
        self.e.setup()
        self.e.get_efficiencies(self.sd, fom_limits, progress_bar=progress_bar)
        self.e.save(detec_tables_dir, model_name)

    @abstractmethod
    def loop(
        self,
        detec_tables_dir: str,
        skip_control_ix: Optional[List] = None,
        **kwargs,
    ):
        """
        Loop over each possible sigma_kern, then each possible peak_appmag, then each row in that corresponding SimDetecTable.
        For each row, inject a Simulation with the specified parameters into a random control light curve.
        Update the row with information about where it was injected, what/where its max FOM is, etc.

        :param detec_tables_dir: Directory where the SimDetecTables should be saved.
        :param skip_control_ix: List of indices of control light curves which may NOT be randomly selected to have a Simulation injected.
        """
        pass


class AtlasSimDetecLoop(SimDetecLoop):
    def __init__(self, sigma_kerns: List, **kwargs):
        super().__init__(sigma_kerns, **kwargs)

    def get_brightness_param_from_detec_tables(
        self, model_name: str, detec_tables_dir: str, param_name="peak_appmag", **kwargs
    ):
        return super().get_brightness_param_from_detec_tables(
            model_name, detec_tables_dir, param_name=param_name, **kwargs
        )

    def get_brightness_param_from_sim_tables(
        self, model_name: str, sim_tables_dir: str, param_name="peak_appmag", **kwargs
    ):
        return super().get_brightness_param_from_sim_tables(
            model_name, sim_tables_dir, param_name=param_name, **kwargs
        )

    def get_max_fom_indices(
        self, sim_lc: SimDetecLightCurve, time_peak_mjd=None, sigma_sim=None, **kwargs
    ):
        """
        Get indices of MJD within 1 sigma of the peak MJD.
        For Gaussians, use the sigma provided.
        For Charlie's model, use the manually calculated value of 2.8.
        """
        if time_peak_mjd is None:
            raise RuntimeError("A peak MJD is required to find the max FOM.")
        if sigma_sim is None:
            # replace with default sigma sim for Charlie's model
            sigma_sim = 2.8

        # measurements within 1 sigma of the peak MJD
        indices = sim_lc.ix_inrange(
            colnames=sim_lc.colnames.mjdbin,
            lowlim=time_peak_mjd - sigma_sim,
            uplim=time_peak_mjd + sigma_sim,
        )
        return indices

    def loop(
        self,
        detec_tables_dir: str,
        skip_control_ix: Optional[List] = None,
        **kwargs,
    ):
        if self.brightness_param is None:
            raise RuntimeError("brightness_param cannot be None")
        if self.sd is None:
            raise RuntimeError("SimDetecTables cannot be None")

        if skip_control_ix is not None and len(skip_control_ix) > 0:
            print(f"\nSkipping control light curve indices: {skip_control_ix}")
            self.sn.remove_lc_indices(skip_control_ix)

        # loop through each rolling sum kernel size
        for sigma_kern in self.sigma_kerns:
            print(
                f"\nUsing rolling sum kernel size sigma_kern={format_float(sigma_kern)} days..."
                "\n-----------------------------------------------------"
            )
            self.sn.apply_rolling_sums(sigma_kern, valid_ix=False, pre_mjd0_ix=False)

            # loop through each possible peak apparent magnitude
            for peak_appmag in self.brightness_param.values:
                sim_detec_table = self.sd.get_table(sigma_kern, peak_appmag)
                sim_detec_table.validate_model_name_col()
                print(
                    f"- Commencing {len(sim_detec_table.t)} simulations for peak brightness of {format_float(peak_appmag)} app mag (= {format_float(mag2flux(peak_appmag))} uJy)..."
                )

                # load the Simulation object based on the data in the first row
                # (we assume here that every row adds the same type of model)
                sim = self.load_sim(dict(sim_detec_table.t.loc[0, :]), verbose=True)

                for i in range(len(sim_detec_table.t)):
                    # pick random control light curve
                    rand_control_index = random.choice(self.sn.control_lc_indices)

                    # add the simulated flux to the chosen control light curve
                    params = sim_detec_table.get_params_at_index(i)
                    sim_lc = self.add_simulation_to_lc(
                        sigma_kern, peak_appmag, rand_control_index, sim, **params
                    )

                    # get the max simulated FOM within certain indices of the light curve
                    indices = self.get_max_fom_indices(sim_lc, **params)
                    max_fom_mjd, max_fom = sim_lc.get_max_fom(indices=indices)

                    # update the corresponding row in the SimDetecTable
                    self.update_sd_row(
                        sigma_kern,
                        peak_appmag,
                        i,
                        rand_control_index,
                        max_fom,
                        max_fom_mjd,
                    )

                self.sd.save_detec_table(sigma_kern, peak_appmag, detec_tables_dir)
                print("\tSuccess")

        if skip_control_ix is not None and len(skip_control_ix) > 0:
            self.sn.add_lc_indices(skip_control_ix)

        print("\nFinished generating all SimDetecTables")


# define command line arguments
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

    return parser


if __name__ == "__main__":
    config = load_config("config.ini")
    args = define_args(config).parse_args()

    print(
        f"\nGenerating SimDetecTables for SN {args.tnsname}, filter {args.filter}, MJD bin size of {format_float(args.mjd_bin_size)} days"
    )
    print(f"Simulations model name: {args.model_name}")
    print(f"Weighted Gaussian rolling sum kernel sizes (days): {args.sigma_kerns}")
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

    simdetec = AtlasSimDetecLoop(args.sigma_kerns)
    simdetec.load_sn(
        config["dir"]["output"],
        colnames,
        args.tnsname,
        args.num_controls + len(args.skip_control_ix),
        mjdbinsize=float(args.mjd_bin_size),
        filt=args.filter,
        flag=hexstring_to_int(config["averaging"]["flag"]),
    )

    sim_tables_dir = get_sim_tables_output_dir(config["dir"]["output"], args.tnsname)
    detec_tables_dir = get_detec_tables_output_dir(
        config["dir"]["output"], args.tnsname
    )

    print()
    simdetec.get_brightness_param_from_sim_tables(args.model_name, sim_tables_dir)
    simdetec.load_sim_tables(args.model_name, sim_tables_dir)
    simdetec.loop(detec_tables_dir, skip_control_ix=args.skip_control_ix)
