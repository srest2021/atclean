#!/usr/bin/env python

"""
Generate a table of simulations for each peak apparent magnitude
using the parameters in simulation_settings.json.
"""

from abc import ABC, abstractmethod
import itertools
import argparse
import os
import sys
import pandas as pd
import numpy as np
from typing import Callable, Dict, List, Optional

from download import make_dir_if_not_exists
from pdastro import pdastrostatsclass
from utils import load_config, load_json_config, abbreviate_list

GAUSSIAN_MODEL_NAME = "gaussian"
ASYMMETRIC_GAUSSIAN_MODEL_NAME = "asymmetric_gaussian"
SIM_TABLE_REQUIRED_COLUMNS = ["model_name", "filename"]


# define command line arguments
def define_args(parser=None, usage=None, conflict_handler="resolve"):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
    parser.add_argument("tnsname", type=str, help="transient name")
    parser.add_argument(
        "-m", "--model_name", type=str, default="gaussian", help="name of model to use"
    )
    parser.add_argument(
        "--sim_config_file",
        default="simulation_settings.json",
        type=str,
        help="file name of JSON file with model information and SimTable generation settings",
    )
    parser.add_argument(
        "--config_file",
        default="config.ini",
        type=str,
        help="file name of .ini file with settings for this class",
    )
    return parser


class Param(ABC):
    """
    Abstract base class for defining simulation parameters.

    Attributes:
        name (str): The name of the parameter.
        values (Optional[List]): The list of values for the parameter.
        is_time_param (bool): Indicates if the parameter is related to time (e.g., MJD).
        is_peak_appmag_param (bool): Indicates if the parameter is the peak apparent magnitude.

    Methods:
        generate(**kwargs): Abstract method to generate parameter values.
        validate_time_param(): Validates and adjusts time parameter values to match the MJDbin format.
        validate_peak_appmag_param(): Validates and adjusts peak apparent magnitude values to two decimal places.
    """

    def __init__(
        self,
        name: str,
        values: Optional[List] = None,
        is_time_param: bool = False,
        is_peak_appmag_param: bool = False,
    ):
        if is_time_param and is_peak_appmag_param:
            raise ValueError(
                "Param cannot be a time parameter and peak apparent magnitude parameter at the same time"
            )

        self.name = name
        self.values = values
        self.is_time_param = is_time_param
        self.is_peak_appmag_param = is_peak_appmag_param

    @abstractmethod
    def generate(self, **kwargs):
        pass

    def validate_time_param(self):
        print(
            f"Making sure the time parameter '{self.name}' values match the MJDbin column format..."
        )
        if self.values:
            self.values = list(np.floor(self.values) + 0.5)

    def validate_peak_appmag_param(self):
        print(
            f"Making sure the peak apparent magnitude parameter '{self.name}' values have up to 2 decimal places..."
        )
        if self.values:
            self.values = [round(v, 2) for v in self.values]

    def __str__(self):
        out = f"Parameter '{self.name}'"
        if self.is_time_param:
            out += " (time param)"
        if self.is_peak_appmag_param:
            out += " (peak apparent magnitude param)"
        out += ": "
        if self.values:
            out += abbreviate_list(self.values)
        else:
            out += "no values yet (call generate() to generate list of values)"
        return out


class ListParam(Param):
    def __init__(
        self,
        name: str,
        values: Optional[List],
        is_time_param: bool = False,
        is_peak_appmag_param: bool = False,
    ):
        super().__init__(
            name,
            values=values,
            is_time_param=is_time_param,
            is_peak_appmag_param=is_peak_appmag_param,
        )

    def generate(self, **kwargs):
        pass


class RangeParam(Param):
    def __init__(
        self,
        name: str,
        minval: float,
        maxval: float,
        step: float,
        is_time_param: bool = False,
        is_peak_appmag_param: bool = False,
    ):
        super().__init__(
            name, is_time_param=is_time_param, is_peak_appmag_param=is_peak_appmag_param
        )
        self.generate(minval, maxval, step)
        if self.is_time_param:
            self.validate_time_param()
        if self.is_peak_appmag_param:
            self.validate_peak_appmag_param()

    def generate(self, minval: float, maxval: float, step: float):
        print(f"Setting to range from {minval} to {maxval} with step size {step}")
        if maxval <= minval:
            raise RuntimeError("Max value must be greater than min value.")
        if step > abs(maxval - minval):
            raise RuntimeError(
                "Step size cannot be greater than the difference between min value and max value."
            )
        self.values = list(np.arange(minval, maxval, step))


class LogRangeParam(Param):
    def __init__(
        self,
        name: str,
        minval: float,
        maxval: float,
        base: int,
        n: int,
        to_int=False,
        is_time_param: bool = False,
        is_peak_appmag_param: bool = False,
    ):
        super().__init__(
            name, is_time_param=is_time_param, is_peak_appmag_param=is_peak_appmag_param
        )
        self.generate(minval, maxval, base, n, to_int=to_int)
        if self.is_time_param:
            self.validate_time_param()
        if self.is_peak_appmag_param:
            self.validate_peak_appmag_param()

    def generate(self, minval: float, maxval: float, base: int, n: int, to_int=False):
        print(
            f'Generating {n}-length {"integer" if to_int else "float"} range using log base {base}'
        )
        if maxval <= minval:
            raise RuntimeError("Max value must be greater than min value.")
        minlog = np.log(minval) / np.log(base)
        maxlog = np.log(maxval) / np.log(base)
        res = list(np.logspace(minlog, maxlog, num=n, base=base))
        if to_int:
            res = [round(num) for num in res]
        self.values = res


class RandomParam(Param):
    def __init__(
        self,
        name: str,
        minval: float,
        maxval: float,
        n: int,
        to_int=False,
        is_time_param: bool = False,
        is_peak_appmag_param: bool = False,
    ):
        super().__init__(
            name, is_time_param=is_time_param, is_peak_appmag_param=is_peak_appmag_param
        )
        self.generate(minval, maxval, n, to_int=to_int)
        if self.is_time_param:
            self.validate_time_param()
        if self.is_peak_appmag_param:
            self.validate_peak_appmag_param()

    def generate(self, minval: float, maxval: float, n: int, to_int=False):
        print(f"Generating {n}-length random list")
        if maxval <= minval:
            raise RuntimeError("maxval must be greater than minval.")
        res = list(np.random.uniform(minval, maxval, n))
        if to_int:
            res = [round(num) for num in res]
        self.values = res


class RandomInRangeParam(Param):
    def __init__(
        self,
        name,
        valid_ranges: List[List[float]],
        n: int,
        is_time_param: bool = False,
        is_peak_appmag_param: bool = False,
    ):
        super().__init__(
            name, is_time_param=is_time_param, is_peak_appmag_param=is_peak_appmag_param
        )
        self.generate(valid_ranges, n)
        if self.is_time_param:
            self.validate_time_param()
        if self.is_peak_appmag_param:
            self.validate_peak_appmag_param()

    def _filter_valid_draws(self, valid_ranges: List[List[float]], draws: List):
        return [
            draw
            for draw in draws
            if any(
                valid_range[0] <= draw <= valid_range[1] for valid_range in valid_ranges
            )
        ]

    def _rec_get_valid_draws(self, valid_ranges: List[List[float]], n: int):
        m = 2 * n
        m_draws = list(np.random.uniform(valid_ranges[0][0], valid_ranges[-1][1], m))
        valid_draws = self._filter_valid_draws(valid_ranges, m_draws)

        m_prime = len(valid_draws)
        if m_prime >= n:
            return valid_draws[:n]
        return valid_draws + self._rec_get_valid_draws(valid_ranges, n - m_prime)

    def generate(self, valid_ranges: List[List[float]], n: int):
        print(
            f"Generating {n}-length random list within the following valid ranges: {valid_ranges}"
        )
        self.values = self._rec_get_valid_draws(valid_ranges, n)


class Params:
    def __init__(self):
        self.d: Dict[str, Param] = {}

    def add(self, param: Param):
        if self.has(param):
            print(f"WARNING: Param {param.name} already exists in list; overwriting...")
        self.d[param.name] = param

    def has(self, param_name):
        return param_name in self.d.keys()

    def validate(self):
        if self.d:
            if not self.has_peak_appmag_param():
                raise RuntimeError(
                    "Peak apparent magnitude parameter missing from parameters"
                )
            if not self.has_time_param():
                raise RuntimeError(f"Time parameter missing from parameters")
        else:
            print("WARNING: Cannot validate params because it is empty")

    def get_num_rows(self):
        total = 1
        for param in self.d.values():
            if param.values and not param.is_peak_appmag_param:
                total *= len(param.values)
        return total

    def has_time_param(self):
        for param in self.d.values():
            if param.is_time_param:
                return True
        return False

    def get_time_param(self):
        for param in self.d.values():
            if param.is_time_param:
                return param
        raise RuntimeError(f"Time parameter missing from parameters: {self.d.keys()}")

    def has_peak_appmag_param(self):
        for param in self.d.values():
            if param.is_peak_appmag_param:
                return True
        return False

    def get_peak_appmag_param(self):
        for param in self.d.values():
            if param.is_peak_appmag_param:
                return param
        raise RuntimeError(
            f"Peak apparent magnitude parameter missing from parameters: {self.d.keys()}"
        )

    def all_names_except_peak_appmag(self):
        return [
            param.name for param in self.d.values() if not param.is_peak_appmag_param
        ]

    def all_names_except_time(self):
        return [param.name for param in self.d.values() if not param.is_time_param]

    def all_params_except_peak_appmag(self):
        return [param for param in self.d.values() if not param.is_peak_appmag_param]

    def all_params_except_time(self):
        return [param for param in self.d.values() if not param.is_time_param]

    def __str__(self):
        out = f"Params list (length {len(self.d)}): "
        for param in self.d.values():
            out += f"\n- {param}"
        return out


def parse_param(
    param_name: str,
    param_info: Dict,
    is_time_param: bool = False,
    is_peak_appmag_param: bool = False,
) -> Param:
    """
    Generate a list of possible values for the parameter using settings from the config file.

    :param_name: Name of parameter as in config file.
    :param_info: Dictionary corresponding to the JSON data under the given parameter in the config file.
    :is_time_param: Is this parameter defining the MJD of the simulation peak, onset, or other time-related property?
    :is_peak_appmag_param: Is this parameter defining the peak apparent magnitude?
    """
    print(f"\nParsing parameter {param_name}:")

    if param_info["type"] == "list":
        res = ListParam(
            param_name,
            param_info["list"],
            is_time_param=is_time_param,
            is_peak_appmag_param=is_peak_appmag_param,
        )

    elif param_info["type"] == "range":
        res = RangeParam(
            param_name,
            param_info["range"]["minval"],
            param_info["range"]["maxval"],
            param_info["range"]["step"],
            is_time_param=is_time_param,
            is_peak_appmag_param=is_peak_appmag_param,
        )

    elif param_info["type"] == "logrange":
        res = LogRangeParam(
            param_name,
            param_info["logrange"]["minval"],
            param_info["logrange"]["maxval"],
            param_info["logrange"]["base"],
            param_info["logrange"]["n"],
            to_int=param_info["logrange"]["to_int"],
            is_time_param=is_time_param,
            is_peak_appmag_param=is_peak_appmag_param,
        )

    elif param_info["type"] == "random":
        res = RandomParam(
            param_name,
            param_info["random"]["minval"],
            param_info["random"]["maxval"],
            param_info["random"]["n"],
            to_int=param_info["random"]["to_int"],
            is_time_param=is_time_param,
            is_peak_appmag_param=is_peak_appmag_param,
        )

    elif param_info["type"] == "random_inrange":
        res = RandomInRangeParam(
            param_name,
            param_info["random_inrange"]["valid_ranges"],
            param_info["random_inrange"]["n"],
            is_time_param=is_time_param,
            is_peak_appmag_param=is_peak_appmag_param,
        )

    else:
        raise RuntimeError(
            "Type must be one of the following: list, range, logrange, random, random_inrange."
        )

    print("Result: ", res.__str__())
    return res


def parse_params(
    model_settings: Dict,
    time_param_name: str = "peak_mjd",
    peak_appmag_param_name: str = "peak_appmag",
) -> Params:
    """
    Parse the parameters in the config file and generate lists of possible values for each parameter.
    """
    params: Params = Params()
    for param_name in model_settings["parameters"]:
        param = parse_param(
            param_name,
            model_settings["parameters"][param_name],
            is_time_param=param_name == time_param_name,
            is_peak_appmag_param=param_name == peak_appmag_param_name,
        )
        params.add(param)
    params.validate()
    print("\n", params)
    return params


class SimTable(pdastrostatsclass):
    def __init__(self, peak_appmag: float, **kwargs):
        """
        Initialize a SimTable.

        :peak_appmag: Peak apparent magnitude for all simulations in this table.
        """
        pdastrostatsclass.__init__(self, **kwargs)
        self.peak_appmag = peak_appmag

    def add_row(self, data: Dict):
        """
        Add a row to the end of the table.

        :param data: Dictionary of column-value pairs.
        """
        self.t = pd.concat([self.t, pd.DataFrame([data])], ignore_index=True)

    def get_sim_filename(self, model_name, tables_dir):
        return f"{tables_dir}/sim_{model_name}_{self.peak_appmag:0.2f}.txt"

    def save_sim_table(self, model_name, tables_dir, verbose=False):
        filename = self.get_sim_filename(model_name, tables_dir)
        if verbose:
            print(f"Saving SimTable {filename}...")
        self.write(filename=filename, overwrite=True, index=False)

    def load_sim_table(self, model_name, tables_dir):
        filename = self.get_sim_filename(model_name, tables_dir)
        try:
            self.load_spacesep(filename, delim_whitespace=True)
        except Exception as e:
            raise RuntimeError(f"Could not load SimTable at {filename}: {str(e)}")

    def __str__(self):
        return self.t.to_string()


class SimTables:
    def __init__(self, model_name: str):
        """
        Initialize a collection of SimTable.

        :param model_name: Name of the model to be used assigned in the config file.
        """
        self.d: Dict[str, SimTable] = {}
        self.peak_appmags = None
        self.model_name = model_name

    def set_peak_appmags(self, peak_appmags: Param):
        if not peak_appmags.is_peak_appmag_param:
            raise ValueError(
                f"peak_appmags.is_peak_appmag_param must be True (got {peak_appmags.is_peak_appmag_param})"
            )

        self.peak_appmags = peak_appmags

    def generate(
        self,
        params: Params,
        filename: str = None,
        mjd_colname: Optional[bool] = False,
        mag_colname: Optional[bool] = False,
        flux_colname: Optional[bool] = False,
    ):
        """
        Generate the table content of each SimTable using parsed parameters from the config file.

        :param params: Params object containing dictionary of Param names and objects.
        :param filename: File name of the model to be used (None if using Gaussian model).
        :param mjd_colname: MJD column name in the model file (None if present but no column name; False if not present).
        :param mag_colname: Magnitude column name in the model file (None if present but no column name; False if not present).
        :param flux_colname: Flux column name in the model file (None if present but no column name; False if not present).
        """
        row = {
            "model_name": self.model_name,
            "filename": np.nan if filename is None else filename,
        }
        if not mjd_colname is False:
            row["mjd_colname"] = np.nan if mjd_colname is None else mjd_colname
        if not mag_colname is False:
            row["mag_colname"] = np.nan if mag_colname is None else mag_colname
        if not flux_colname is False:
            row["flux_colname"] = np.nan if flux_colname is None else flux_colname

        print()
        self.d = {}
        num_rows = params.get_num_rows()
        self.set_peak_appmags(params.get_peak_appmag_param())
        for peak_appmag in self.peak_appmags.values:
            print(
                f"Generating {num_rows}-length SimTable for peak_appmag={peak_appmag}..."
            )
            self.d[peak_appmag] = SimTable(peak_appmag)

            combinations = list(
                itertools.product(
                    *(param.values for param in params.all_params_except_peak_appmag())
                )
            )
            combinations_dicts = [
                dict(zip(params.all_names_except_peak_appmag(), combo))
                for combo in combinations
            ]

            for combination in combinations_dicts:
                combination.update({"peak_appmag": peak_appmag})
                combination.update(row)
                self.d[peak_appmag].add_row(combination)

        print("Success")

    def save_all(self, tables_dir: str):
        print(f"\nSaving SimTables in directory: {tables_dir}")

        if self.peak_appmags is None:
            raise RuntimeError(
                "Cannot save SimTables: missing peak apparent magnitudes"
            )

        if not self.d:
            raise RuntimeError(
                "Cannot save SimTables: empty dictionary (call generate() first)"
            )

        make_dir_if_not_exists(tables_dir)

        for peak_appmag in self.peak_appmags.values:
            self.d[peak_appmag].save_sim_table(self.model_name, tables_dir)
        print("Success")

    def load_all(self, tables_dir: str, peak_appmags: Param):
        print(f"\nLoading SimTables in directory: {tables_dir}")
        self.d = {}
        self.set_peak_appmags(peak_appmags)

        for peak_appmag in self.peak_appmags.values:
            self.d[peak_appmag] = SimTable(peak_appmag)
            self.d[peak_appmag].load_sim_table(self.model_name, tables_dir)
        print("Success")


def get_sim_tables_output_dir(output_dir: str, tnsname: str):
    return os.path.join(output_dir, tnsname, "bump_analysis", "sim_tables")


def parse_colname_info(model_settings: Dict, model_name: str):
    filename, mjd_colname, mag_colname, flux_colname = None, False, False, False
    if (
        model_name != GAUSSIAN_MODEL_NAME
        and model_name != ASYMMETRIC_GAUSSIAN_MODEL_NAME
    ):
        try:
            filename = model_settings["filename"]
            mjd_colname = model_settings["mjd_column_name"]
            mag_colname = model_settings["mag_column_name"]
            flux_colname = model_settings["flux_column_name"]

            if " " in filename:
                raise RuntimeError("Filename cannot have spaces.")

            if mag_colname is False and flux_colname is False:
                raise RuntimeError(
                    f"Model must have either mag or flux column. Please set one or both fields to null or the correct column name."
                )
        except Exception as e:
            raise RuntimeError(f"{str(e)}")

    return filename, mjd_colname, mag_colname, flux_colname


if __name__ == "__main__":
    args = define_args().parse_args()
    config = load_config(args.config_file)
    sim_config = load_json_config(args.sim_config_file)

    if " " in args.model_name:
        raise RuntimeError("Model name cannot have spaces.")
    if args.model_name not in sim_config:
        raise RuntimeError(
            f"Model '{args.model_name}' not found in simulation config file\n"
            f"Available models: {', '.join(sim_config.keys())}"
        )
    print(f"Loading settings for model '{args.model_name}'...")
    model_settings = sim_config[args.model_name]

    print("Parsing model parameters...")
    params = parse_params(
        model_settings,
        time_param_name=model_settings["time_parameter_name"],
        peak_appmag_param_name=model_settings["peak_appmag_parameter_name"],
    )
    filename, mjd_colname, mag_colname, flux_colname = parse_colname_info(
        model_settings, args.model_name
    )

    sim_tables = SimTables(args.model_name)
    sim_tables.generate(
        params,
        filename=filename,
        mjd_colname=mjd_colname,
        mag_colname=mag_colname,
        flux_colname=flux_colname,
    )
    sim_tables_output_dir = get_sim_tables_output_dir(
        config["dir"]["output"], args.tnsname
    )
    sim_tables.save_all(sim_tables_output_dir)
