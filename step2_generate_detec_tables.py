#!/usr/bin/env python

from abc import ABC, abstractmethod
from collections import defaultdict
from configparser import ConfigParser
import itertools
import json
import os
import random
import argparse, re
from copy import deepcopy
import sys
from typing import Dict, List, Optional, Self, Tuple, Union
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
    get_sim_tables_output_dir,
    GAUSSIAN_MODEL_NAME,
    ASYMMETRIC_GAUSSIAN_MODEL_NAME,
)
from lightcurve import SimDetecLightCurve, SimDetecSupernova, Simulation
from utils import (
    PresetColumnNames,
    extract_from_subdir,
    format_float,
    get_allowed_presets,
    hexstring_to_int,
    load_config,
)


NON_PARAM_COLNAMES = {
    "sigma_kern",
    "filename",
    "model_name",
    "mjd_colname",
    "mag_colname",
    "flux_colname",
    "control_index",
    "filter",
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


def get_brightness_values_from_dir(directory: str, pattern: re.Pattern):
    brightnesses = extract_from_subdir(directory, pattern, 1, convert_function=float)
    res = list(brightnesses)
    res.sort()
    return res


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

    def get_param_colnames(
        self, skip_time_col: bool = False, skip_brightness_col: bool = False
    ) -> List[str]:
        return [
            col
            for col in self.t.columns
            if not (
                (skip_time_col and col.startswith(TIME_PARAM_PREFIX))
                or (skip_brightness_col and col.startswith(BRIGHTNESS_PARAM_PREFIX))
                or col in NON_PARAM_COLNAMES
            )
        ]

    def get_params_at_index(
        self, index: int, skip_time_col: bool = False, skip_brightness_col: bool = True
    ) -> Dict:
        """
        Get a dictionary of the parameter column-value pairs of the Simulation object at a certain row.
        Any known non-parameter column names (including brightness and, optionally, time) will be skipped.

        :param index: Index of the table from which to get the parameter column-value pairs.
        :param skip_time_col: Whether to exclude time columns (i.e., those starting with TIME_PARAM_PREFIX).
        """
        colnames = self.get_param_colnames(
            skip_time_col=skip_time_col, skip_brightness_col=skip_brightness_col
        )
        return dict(self.t.loc[index, colnames])

    def get_params(
        self, skip_time_col: bool = True, skip_brightness_col: bool = True
    ) -> Params:
        params = Params()

        colnames = self.get_param_colnames(
            skip_time_col=skip_time_col, skip_brightness_col=skip_brightness_col
        )
        for col in colnames:
            values = self.t[col].unique()
            param = ListParam(col, values)
            params.add(param)

        return params

    def update_row(
        self,
        index: int,
        filt: str,
        control_index: int,
        max_fom: float,
        max_fom_mjd: float,
    ):
        """
        Update a certain row of the table.

        :param index: Index of the row to update.
        :param filt: Filter of the light curve into which the simulation was injected.
        :param control_index: Control index of the light curve into which the simulation was injected.
        :param max_fom: Maximum FOM value of the simulated light curve within a certain range of the injection time.
        :param max_fom_mjd: MJD of the max_fom value.
        """
        data = {
            "control_index": control_index,
            "filter": filt,
            "max_fom": max_fom,
            "max_fom_mjd": max_fom_mjd,
        }
        for key, value in data.items():
            self.t.at[index, key] = value

    def get_detec_filename(
        self, model_name: str, filt: str, detec_tables_dir: str
    ) -> str:
        """
        Get the filename of the SimDetecTable.

        :param model_name: Name of the model of which the SimDetecTable contains simulations.
        :param detec_tables_dir: Directory where the SimDetecTable is located.
        """
        return f"{detec_tables_dir}/simdetec_{model_name}_{format_float(self.sigma_kern)}_{format_float(self.brightness)}.{filt}.txt"

    def load_detec_table(self, model_name: str, filt: str, detec_tables_dir: str):
        """
        Load an existing SimDetecTable.

        :param model_name: Name of the model of which the SimDetecTable contains simulations.
        :param detec_tables_dir: Directory where the SimDetecTable is located.
        """
        filename = self.get_detec_filename(model_name, filt, detec_tables_dir)
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

    def save_detec_table(self, model_name: str, filt: str, detec_tables_dir: str):
        """
        Save the current SimDetecTable.

        :param model_name: Name of the model of which the SimDetecTable contains simulations.
        :param detec_tables_dir: Directory where the SimDetecTable should be saved.
        """
        filename = self.get_detec_filename(model_name, filt, detec_tables_dir)
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
        self,
        filt: str,
        brightness_param: Param,
        model_name: str,
        sigma_kerns: List[float],
    ):
        self.filt: str = filt
        self.model_name: str = model_name
        self.sigma_kerns: List[float] = sigma_kerns
        self.brightness_param = brightness_param
        self.d: Dict[float, Dict[float, SimDetecTable]] = {}

    def get_table(self, sigma_kern: float, brightness: float):
        return self.d[sigma_kern][brightness]

    def update_row(
        self,
        sigma_kern: float,
        brightness: float,
        index: int,
        filt: str,
        control_index: int,
        max_fom: float,
        max_fom_mjd: float,
    ):
        """
        Update a certain row of a SimDetecTable.

        :param sigma_kern: Rolling sum kernel size corresponding to the SimDetecTable to update.
        :param brightness: Brightness corresponding to the SimDetecTable to update.
        :param index: Index of the row to update in the SimDetecTable.
        :param filt: Filter of the light curve into which the simulation was injected.
        :param control_index: Control index of the light curve into which the simulation was injected.
        :param max_fom: Maximum FOM value of the simulated light curve within a certain range of the injection time.
        :param max_fom_mjd: MJD of the max_fom value.
        """
        self.d[sigma_kern][brightness].update_row(
            index, filt, control_index, max_fom, max_fom_mjd
        )

    def get_efficiency(
        self, sigma_kern: float, brightness: float, fom_limit: float, **params
    ):
        return self.d[sigma_kern][brightness].get_efficiency(fom_limit, **params)

    def save_detec_table(
        self, sigma_kern: float, brightness: float, detec_tables_dir: str
    ):
        self.d[sigma_kern][brightness].save_detec_table(
            self.model_name, self.filt, detec_tables_dir
        )

    def save_all(self, detec_tables_dir: str):
        print(f"\nSaving SimDetecTables in directory: {detec_tables_dir}")
        make_dir_if_not_exists(detec_tables_dir)
        for sigma_kern in self.d.keys():
            for table in self.d[sigma_kern].values():
                table.save_detec_table(self.model_name, self.filt, detec_tables_dir)
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
                    self.model_name, self.filt, detec_tables_dir
                )
        print("Success")


class SimulationFactory:
    def __init__(self, verbose: bool = False):
        self.verbose = verbose

    def _log(self, msg: str):
        if self.verbose:
            print(msg)

    def _parse_colname_val(self, colname: str, row: dict):
        if colname in row:
            return None if np.isnan(row[colname]) else row[colname]
        return False

    def from_table(
        self, table: Union[SimTable, SimDetecTable], row_index: int
    ) -> Simulation:
        """
        Create a Simulation object given a row from a SimTable or SimDetecTable.
        """
        row = dict(table.t.loc[row_index, :])
        model_name = row["model_name"]
        filename = None if not isinstance(row["filename"], str) else row["filename"]

        mjd_colname = self._parse_colname_val("mjd_colname", row)
        mag_colname = self._parse_colname_val("mag_colname", row)
        flux_colname = self._parse_colname_val("flux_colname", row)

        if model_name == GAUSSIAN_MODEL_NAME:
            self._log("\tConstructing Gaussian simulation")
            return Gaussian()
        elif model_name == ASYMMETRIC_GAUSSIAN_MODEL_NAME:
            self._log("\tConstructing AsymmetricGaussian simulation")
            return AsymmetricGaussian()
        else:
            self._log(
                f"\tConstructing '{model_name}' simulation with MJD column {mjd_colname}, mag column {mag_colname}, flux column {flux_colname}, filename: {filename}"
            )
            return Model(
                filename=filename,
                mjd_colname=mjd_colname,
                mag_colname=mag_colname,
                flux_colname=flux_colname,
                model_name=model_name,
            )

    @staticmethod
    def new_gaussian() -> Simulation:
        return Gaussian()

    @staticmethod
    def new_asymmetric_gaussian() -> Simulation:
        return AsymmetricGaussian()


class InjectionLoop(ABC):
    def __init__(
        self,
        sigma_kerns: List[float],
        model_name: str,
        sim_tables_dir: str,
        detec_tables_dir: str,
        **kwargs,
    ):
        self.sigma_kerns: List[float] = sigma_kerns
        self.model_name = model_name
        self.sim_tables_dir = sim_tables_dir
        self.detec_tables_dir = detec_tables_dir

        self.brightness_param: Param = None
        self.sn: SimDetecSupernova = None
        self.tables: SimDetecTables = None

    def get_brightness_param_from_sim_tables(
        self,
        param_name: str = "brightness",
    ):
        pattern = re.compile(rf"^sim_{re.escape(self.model_name)}_(\d+\.\d+)\.txt$")
        values = get_brightness_values_from_dir(self.sim_tables_dir, pattern)
        self.brightness_param = ListParam(
            param_name, values, param_type=ParamType.BRIGHTNESS
        )
        print(self.brightness_param)

    def set_brightness_param(self, values: List[float], param_name="brightness"):
        self.brightness_param = ListParam(
            param_name, values, param_type=ParamType.BRIGHTNESS
        )

    def load_sn(
        self,
        data_dir: str,
        colnames: PresetColumnNames,
        tnsname: str,
        num_controls: int,
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
        :param mjdbinsize: MJD bin size of the averaged light curves to load.
        :param filt: Filter of the averaged light curves to load.
        :param mjd_ranges: Valid MJD ranges into which we inject Simulations.
        :param skip_control_ix: List of indices of control light curves which may NOT be randomly selected to have a Simulation injected.
        :param flag: Flag that denotes bad days in the binned light curves to load.
        """
        self.sn = SimDetecSupernova(
            colnames, tnsname, mjdbinsize=mjdbinsize, filt=filt, flag=flag
        )
        self.sn.load_all(data_dir, num_controls=num_controls)
        self.sn.remove_rolling_sums()
        self.sn.remove_simulations()
        if mjd_ranges is not None:
            self.sn.set_mjd_ranges(mjd_ranges)
        if skip_control_ix is not None and len(skip_control_ix) > 0:
            print(f"Skipping control light curve indices: {skip_control_ix}")
            self.sn.remove_lc_indices(skip_control_ix)

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

        self.sn = sn
        self.sn.remove_rolling_sums()
        self.sn.remove_simulations()
        if mjd_ranges is not None:
            self.sn.set_mjd_ranges(mjd_ranges)
        if skip_control_ix is not None and len(skip_control_ix) > 0:
            print(f"Skipping control light curve indices: {skip_control_ix}")
            self.sn.remove_lc_indices(skip_control_ix)

    def load_sim_tables(self):
        """
        Load existing SimTables and construct SimDetecTables out of them.
        """
        self.tables = SimDetecTables(
            self.sn.filt, self.brightness_param, self.model_name, self.sigma_kerns
        )
        self.tables.load_all_from_sim_tables(self.sim_tables_dir)

    def load_detec_tables(self):
        """
        Load existing SimDetecTables.
        """
        self.tables = SimDetecTables(
            self.sn.filt, self.brightness_param, self.model_name, self.sigma_kerns
        )
        self.tables.load_all(self.detec_tables_dir)

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
        :param control_index: The control index of the light curve to add the Simulation to.
        :param sim: The Simulation to add.
        :param remove_old: Remove any old simulations before adding the new simulated flux.
        :param kwargs: Additional Simulation parameters (e.g., sigma_sim=1.0 and time_peak_mjd=56780.5 for Gaussian)
        """
        if verbose:
            print(f"Adding simulation: {sim}")
            if kwargs:
                print(f"Additional simulation parameters: {kwargs}")
        if self.sn is None:
            raise ValueError(
                "SN must be set via self.load_sn() or self.set_sn() before injecting a simulation"
            )
        if control_index not in self.sn.control_lc_indices:
            raise ValueError(
                f"Control index {control_index} missing from control light curve indices (may have been excluded): {self.sn.control_lc_indices}"
            )

        lc = deepcopy(self.sn.lcs[control_index])
        if not lc.colnames.snrsumnorm in lc.t.columns:
            print(
                "WARNING: Rolling sum not applied to light curve prior to injecting simulation"
            )

        good_ix = lc.get_good_indices(flag=self.sn.flag)

        sim_flux = sim.get_sim_flux(
            lc.t.loc[good_ix, self.sn.colnames.mjdbin], brightness, **kwargs
        )

        lc.add_sim_flux(
            good_ix,
            sim_flux,
            cur_sigma_kern=sigma_kern,
            verbose=verbose,
            remove_old=remove_old,
        )
        return sim_flux, lc

    @abstractmethod
    def get_injection_search_indices(
        self,
        sim_lc: SimDetecLightCurve,
        **kwargs,
    ):
        """
        From a light curve with a Simulation injected, return the indices within which to search for the max FOM.

        :param sim_lc: SimDetecLightCurve with a Simulation injected.
        :param kwargs: Any additional helpful parameters (e.g., peak MJD and the sigma or size of the Simulation) that can be used to calculate indices within a certain range of the injection time.
        """
        pass

    def inject_and_record_sim(
        self,
        sim_detec_table: SimDetecTable,
        row_index: int,
        sim: Simulation,
        sigma_kern: float,
        brightness: float,
    ):
        # pick random control light curve
        rand_control_index = random.choice(self.sn.control_lc_indices)

        # add the simulated flux to the chosen control light curve
        params = sim_detec_table.get_params_at_index(row_index)
        _, sim_lc = self.add_simulation_to_lc(
            sigma_kern, brightness, rand_control_index, sim, **params
        )

        # get the max simulated FOM within certain indices of the light curve
        indices = self.get_injection_search_indices(sim_lc, **params)
        max_fom_mjd, max_fom = sim_lc.get_max_fom(indices=indices)

        # update the corresponding row in the SimDetecTable
        self.tables.update_row(
            sigma_kern,
            brightness,
            row_index,
            self.sn.filt,
            rand_control_index,
            max_fom,
            max_fom_mjd,
        )

    # def calculate_efficiencies(
    #     self,
    #     fom_limits: (
    #         List[float]
    #         | List[List[float]]
    #         | Dict[float, float]
    #         | Dict[float, List[float]]
    #     ),
    #     params: Params,
    #     detec_tables_dir: str,
    #     model_name: str,
    #     progress_bar: bool = True,
    # ):
    #     """
    #     Construct and save an EfficiencyTable that contains efficiencies for every combination of a Simulation's sigma_kern, peak_appmag, and other parameters EXCEPT the time parameter.

    #     :param fom_limits: Dict or List of FOM limits.
    #     :param params: Collection of parameter names and possible values.
    #     :param detec_tables_dir: Directory where the EfficiencyTable should be saved.
    #     :param model_name: Name of the model for which to calculate efficiencies.
    #     """
    #     self.e = EfficiencyTable(self.sigma_kerns, params)
    #     self.e.setup()
    #     self.e.get_efficiencies(self.sd, fom_limits, progress_bar=progress_bar)
    #     self.e.save(detec_tables_dir, model_name)

    @abstractmethod
    def loop(
        self,
        **kwargs,
    ):
        """
        Loop over each possible sigma_kern, then each possible peak_appmag, then each row in that corresponding SimDetecTable.
        For each row, inject a Simulation with the specified parameters into a random control light curve.
        Update the row with information about where it was injected, what/where its max FOM is, etc.
        """
        pass


class AtlasInjectionLoop(InjectionLoop):
    def __init__(self, sigma_kerns: List, **kwargs):
        super().__init__(sigma_kerns, **kwargs)

    def get_brightness_param_from_sim_tables(self):
        return super().get_brightness_param_from_sim_tables(param_name="peak_appmag")

    def get_injection_search_indices(
        self, sim_lc: SimDetecLightCurve, time_peak_mjd=None, sigma_sim=None
    ):
        """
        Get indices of measurements within 1 sigma of the peak MJD.
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

    def loop(self):
        if self.brightness_param is None:
            raise RuntimeError(
                "self.brightness_param must be set before calling self.loop()"
            )
        if self.tables is None:
            raise RuntimeError(
                "SimDetecTables (self.tables) must be set before calling self.loop()"
            )
        if self.sn is None:
            raise RuntimeError(
                "Supernova (self.sn) must be set before calling self.loop()"
            )

        # loop through each rolling sum kernel size
        for sigma_kern in self.sigma_kerns:
            print(
                f"\nUsing rolling sum kernel size sigma_kern={format_float(sigma_kern)} days..."
                "\n-----------------------------------------------------"
            )
            self.sn.apply_rolling_sums(sigma_kern, valid_ix=True, pre_mjd0_ix=False)

            sim_factory = SimulationFactory(verbose=True)

            # loop through each possible peak apparent magnitude
            for peak_appmag in self.brightness_param.values:
                sim_detec_table = self.tables.get_table(sigma_kern, peak_appmag)
                sim_detec_table.validate_model_name_col()
                print(
                    f"- Commencing {len(sim_detec_table.t)} simulations for peak brightness of {format_float(peak_appmag)} app mag (= {format_float(mag2flux(peak_appmag))} uJy)..."
                )

                # load the Simulation object based on the data in the first row
                # (we assume here that every row adds the same type of model)
                sim = sim_factory.from_table(sim_detec_table, 0)

                for i in range(len(sim_detec_table.t)):
                    self.inject_and_record_sim(
                        sim_detec_table, i, sim, sigma_kern, peak_appmag
                    )

                self.tables.save_detec_table(
                    sigma_kern, peak_appmag, self.detec_tables_dir
                )
                print("\tSuccess")
        print("\nFinished generating all SimDetecTables")


def mjd_range_type(value):
    try:
        ranges = json.loads(value)
    except json.JSONDecodeError:
        raise argparse.ArgumentTypeError("MJD ranges must be valid JSON.")

    if not isinstance(ranges, list):
        raise argparse.ArgumentTypeError("MJD ranges must be a list of lists.")

    for i, r in enumerate(ranges):
        if not (isinstance(r, list) and len(r) == 2):
            raise argparse.ArgumentTypeError(
                f"Range {i} must be a 2-element list, got: {r}"
            )
        if not all(isinstance(x, (int, float, np.integer, np.floating)) for x in r):
            raise argparse.ArgumentTypeError(
                f"Range {i} must contain only numbers, got: {r}"
            )
        if r[0] > r[1]:
            raise argparse.ArgumentTypeError(f"Range {i} has start > end: {r}")
    return ranges


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

    sim_tables_dir = get_sim_tables_output_dir(config["dir"]["output"], args.tnsname)
    detec_tables_dir = get_detec_tables_output_dir(
        config["dir"]["output"], args.tnsname
    )
    injection_loop = AtlasInjectionLoop(
        args.sigma_kerns, args.model_name, sim_tables_dir, detec_tables_dir
    )

    injection_loop.load_sn(
        config["dir"]["output"],
        colnames,
        args.tnsname,
        args.num_controls + len(args.skip_control_ix),
        mjdbinsize=float(args.mjd_bin_size),
        filt=args.filter,
        mjd_ranges=args.mjd_ranges,
        skip_control_ix=args.skip_control_ix,
        flag=hexstring_to_int(config["averaging"]["flag"]),
    )
    if args.mjd_ranges is not None:
        print(f"\nValid MJD ranges: {args.mjd_ranges}")

    print()
    injection_loop.get_brightness_param_from_sim_tables()
    injection_loop.load_sim_tables()
    injection_loop.loop()
