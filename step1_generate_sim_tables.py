#!/usr/bin/env python

"""
Generate a table of simulations for each brightness
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
from enum import Enum, auto

from download import make_dir_if_not_exists
from pdastro import pdastrostatsclass
from utils import (
    CustomLogger,
    format_float_string,
    load_config,
    load_json_config,
    abbreviate_list,
)

GAUSSIAN_MODEL_NAME = "gaussian"
ASYMMETRIC_GAUSSIAN_MODEL_NAME = "asymmetric_gaussian"
TIME_PARAM_PREFIX = "time_"
BRIGHTNESS_PARAM_PREFIX = "brightness_"


# define command line arguments
def define_args(parser=None, usage=None, conflict_handler="resolve"):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
    parser.add_argument("tnsname", type=str, help="transient name")
    parser.add_argument(
        "model_name", type=str, default="gaussian", help="name of model to use"
    )
    parser.add_argument(
        "--step1_config_file",
        default="step1_settings.json",
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


def find_prefix_in_list(l: List[str], prefix: str):
    for item in l:
        if item.startswith(prefix):
            return item
    return None


def remove_prefix(string: str, prefix: str):
    if string.startswith(prefix):
        return string[len(prefix) :]
    return string


def remove_any_prefix(string: str):
    if string.startswith(BRIGHTNESS_PARAM_PREFIX):
        return remove_prefix(string, BRIGHTNESS_PARAM_PREFIX)
    if string.startswith(TIME_PARAM_PREFIX):
        return remove_prefix(string, TIME_PARAM_PREFIX)
    return string


class ParamType(Enum):
    # Indicates if the parameter is related to time (e.g., peak or onset MJD).
    TIME = auto()

    # Indicates if the parameter is related to brightness (e.g., peak magnitude or flux).
    BRIGHTNESS = auto()

    OTHER = auto()


class Param(ABC):
    """
    Abstract base class for defining simulation parameters.
    """

    def __init__(
        self,
        name: str,
        values: Optional[List] = None,
        param_type: ParamType = ParamType.OTHER,
        verbose: bool = True,
    ):
        """
        :param name (str): The name of the parameter.
        :param values (Optional[List]): The list of values for the parameter.
        :param param_type (ParamType): ParamType indicating whether the Param is related to time, brightness, or neither.
        """
        self.logger = CustomLogger(self.__class__.__name__)

        self.param_type: ParamType = param_type
        self.verbose = verbose

        out = f"Creating parameter '{name}'"
        if self.is_time_param:
            out += " (time param)"
        if self.is_brightness_param:
            out += " (brightness param)"
        self._log(out)

        self.name = name
        self._values: Optional[List] = list(values) if values is not None else values
        self.validate_name()

    def _log(self, message: str):
        if self.verbose:
            self.logger.body(message)

    @property
    def values(self) -> List:
        if self._values is None:
            raise RuntimeError("Parameter values cannot be None")
        else:
            return self._values

    @property
    def is_time_param(self) -> bool:
        return self.param_type == ParamType.TIME

    @property
    def is_brightness_param(self) -> bool:
        return self.param_type == ParamType.BRIGHTNESS

    @abstractmethod
    def generate(self, **kwargs):
        """
        Abstract method to generate parameter values.
        """
        pass

    def trim(
        self, min_value: Optional[float] = None, max_value: Optional[float] = None
    ):
        if self._values is None:
            return

        self._values = [
            v
            for v in self._values
            if (min_value is None or v >= min_value)
            and (max_value is None or v <= max_value)
        ]

    def round_to(self, n_digits: int):
        if self._values:
            self._values = [round(v, n_digits) for v in self.values]

    def has_values(self) -> bool:
        return self._values is not None and len(self._values) > 0

    def validate_time_param(self):
        """
        Validates and adjusts time parameter values to match the MJDbin format.
        """
        self._log(
            f"Making sure the time parameter '{self.name}' values match the MJDbin column format"
        )
        if self._values:
            self._values = list(np.floor(self.values) + 0.5)

    def validate_brightness_param(self):
        """
        Validates and adjusts brightness values to two decimal places.
        """
        self._log(
            f"Making sure the brightness parameter '{self.name}' values have up to 2 decimal places"
        )
        if self._values:
            self._values = [round(v, 2) for v in self.values]

    def validate_name(self):
        if self.is_time_param and not self.name.startswith(TIME_PARAM_PREFIX):
            self.name = TIME_PARAM_PREFIX + self.name
            return

        if self.is_brightness_param and not self.name.startswith(
            BRIGHTNESS_PARAM_PREFIX
        ):
            self.name = BRIGHTNESS_PARAM_PREFIX + self.name
            return

    def __str__(self):
        out = f"Parameter '{self.name}'"
        if self.is_time_param:
            out += " (time param)"
        if self.is_brightness_param:
            out += " (brightness param)"
        out += ": "
        if self._values is not None:
            out += abbreviate_list(self.values)
        else:
            out += "no values yet (call generate() to generate list of values)"
        return out

    def __eq__(self, other):
        if not isinstance(other, Param):
            return False

        if self._values is None and other._values is None:
            values_equal = True
        elif self._values is None or other._values is None:
            return False
        else:
            values_equal = np.array_equal(self._values, other._values)

        return (
            self.name == other.name
            and self.param_type == other.param_type
            and values_equal
        )


class ListParam(Param):
    """
    A parameter defined by a fixed list of values.
    """

    def __init__(
        self,
        name: str,
        values: Optional[List],
        param_type: ParamType = ParamType.OTHER,
        verbose: bool = True,
    ):
        """
        :param name (str): The name of the parameter.
        :param values (Optional[List]): The list of values for the parameter.
        :param param_type (ParamType): ParamType indicating whether the Param is related to time, brightness, or neither.
        """
        super().__init__(name, values=values, param_type=param_type, verbose=verbose)

    def generate(self, **kwargs):
        pass


class RangeParam(Param):
    """
    A parameter defined by a range of values with a fixed step size.
    """

    def __init__(
        self,
        name: str,
        minval: float,
        maxval: float,
        step: float,
        param_type: ParamType = ParamType.OTHER,
        verbose: bool = True,
    ):
        """
        :param name (str): The name of the parameter.
        :param minval (float): The minimum value of the range.
        :param maxval (float): The maximum value of the range.
        :param step (float): The step size between consecutive values in the range.
        :param param_type (ParamType): ParamType indicating whether the Param is related to time, brightness, or neither.
        """
        super().__init__(name, param_type=param_type, verbose=verbose)
        self.generate(minval, maxval, step)
        if self.is_time_param:
            self.validate_time_param()
        if self.is_brightness_param:
            self.validate_brightness_param()

    def generate(self, minval: float, maxval: float, step: float):
        self._log(f"Setting to range from {minval} to {maxval} with step size {step}")
        if maxval <= minval:
            raise RuntimeError("Max value must be greater than min value.")
        if step > abs(maxval - minval):
            raise RuntimeError(
                "Step size cannot be greater than the difference between min value and max value."
            )
        self._values = list(np.arange(minval, maxval + step, step))


class LogRangeParam(Param):
    """
    A parameter defined by a logarithmic range of values.
    """

    def __init__(
        self,
        name: str,
        minval: float,
        maxval: float,
        base: int,
        n: int,
        n_digits: int = 5,
        param_type: ParamType = ParamType.OTHER,
        verbose: bool = True,
    ):
        """
        :param name (str): The name of the parameter.
        :param minval (float): The minimum value of the range.
        :param maxval (float): The maximum value of the range.
        :param base (int): The logarithmic base to use.
        :param n (int): The number of values to generate in the range.
        :param n_digits (int): The number of decimal places to round to.
        :param param_type (ParamType): ParamType indicating whether the Param is related to time, brightness, or neither.
        """
        super().__init__(name, param_type=param_type, verbose=verbose)
        self.generate(minval, maxval, base, n, n_digits=n_digits)
        if self.is_time_param:
            self.validate_time_param()
        if self.is_brightness_param:
            self.validate_brightness_param()

    def generate(
        self,
        minval: float,
        maxval: float,
        base: int,
        n: int,
        n_digits: int = 5,
    ):
        self._log(
            f"Generating {n}-length log range of floats rounded to {n_digits} decimal places using log base {base}"
        )
        if maxval <= minval:
            raise RuntimeError("Max value must be greater than min value.")
        minlog = np.log(minval) / np.log(base)
        maxlog = np.log(maxval) / np.log(base)
        res = list(np.logspace(minlog, maxlog, num=n, base=base))
        self._values = res
        self.round_to(n_digits)


class RandomParam(Param):
    """
    A parameter defined by a random set of values within a specified range.
    """

    def __init__(
        self,
        name: str,
        minval: float,
        maxval: float,
        n: int,
        n_digits: int = 5,
        param_type: ParamType = ParamType.OTHER,
        verbose: bool = True,
    ):
        """
        :param name (str): The name of the parameter.
        :param minval (float): The minimum value of the range.
        :param maxval (float): The maximum value of the range.
        :param n (int): The number of random values to generate.
        :param n_digits (int): The number of decimal places to round to.
        :param param_type (ParamType): ParamType indicating whether the Param is related to time, brightness, or neither.
        """
        super().__init__(name, param_type=param_type, verbose=verbose)
        self.generate(minval, maxval, n, n_digits=n_digits)
        if self.is_time_param:
            self.validate_time_param()
        if self.is_brightness_param:
            self.validate_brightness_param()

    def generate(
        self,
        minval: float,
        maxval: float,
        n: int,
        n_digits: int = 5,
    ):
        self._log(
            f"Generating {n}-length random list of floats rounded to {n_digits} decimal places"
        )
        if maxval <= minval:
            raise RuntimeError("maxval must be greater than minval.")
        res = list(np.random.uniform(minval, maxval, n))
        self._values = res
        self.round_to(n_digits)


class RandomInRangeParam(Param):
    """
    A parameter defined by a random set of values within specified valid ranges.
    """

    def __init__(
        self,
        name,
        valid_ranges: List[List[float]],
        n: int,
        param_type: ParamType = ParamType.OTHER,
        verbose: bool = True,
    ):
        """
        :param name (str): The name of the parameter.
        :param valid_ranges (List[List[float]]): A list of valid ranges, where each range is a list of two floats [min, max].
        :param n (int): The number of random values to generate.
        :param param_type (ParamType): ParamType indicating whether the Param is related to time, brightness, or neither.
        """
        super().__init__(name, param_type=param_type, verbose=verbose)
        self.generate(valid_ranges, n)
        if self.is_time_param:
            self.validate_time_param()
        if self.is_brightness_param:
            self.validate_brightness_param()

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
        self._log(
            f"Generating {n}-length random list of floats within the following valid ranges: {valid_ranges}"
        )
        self._values = self._rec_get_valid_draws(valid_ranges, n)


class Params:
    """
    A collection of simulation parameters.
    """

    def __init__(self):
        self.logger = CustomLogger(self.__class__.__name__)
        self.time_param: Optional[Param] = None
        self.brightness_param: Optional[Param] = None
        self.other: Dict[str, Param] = {}

    def add(self, param: Param):
        """
        Adds a parameter to the collection.
        """
        param.validate_name()

        if param.is_time_param:
            if self.time_param is not None:
                self.logger.warning("Time parameter already set")
            self.time_param = param

        elif param.is_brightness_param:
            if self.brightness_param is not None:
                self.logger.warning("Brightness parameter already set")
            self.brightness_param = param

        else:
            if self.has(param.name):
                self.logger.warning(
                    f"Param {param.name} already exists in list; overwriting", dots=True
                )
            self.other[param.name] = param

    def has(self, param_name):
        """
        Checks if a parameter exists in the collection.
        """
        return (
            (self.time_param and self.time_param.name == param_name)
            or (self.brightness_param and self.brightness_param.name == param_name)
            or param_name in self.other.keys()
        )

    def validate(self, check_brightness: bool = True, check_time: bool = True):
        """
        Validates the collection, ensuring required parameters (time and brightness) are present.
        """
        if check_brightness and not self.has_brightness_param():
            raise RuntimeError("Brightness parameter missing from parameters")

        if check_time and not self.has_time_param():
            raise RuntimeError(f"Time parameter missing from parameters")

    def get_num_combinations(self, except_brightness: bool = True) -> int:
        """
        Calculates the total number of rows in the simulation table based on parameter combinations.
        """
        total = 1
        params = (
            self.all_params_except_brightness()
            if except_brightness
            else self.all_params()
        )
        for param in params:
            if param.has_values():
                total *= len(param.values)
        return total

    def get_combinations(self, except_brightness: bool = True) -> List:
        """
        Generates all possible combinations of parameter values, excluding the brightness parameter.
        This method uses the Cartesian product of the values of all parameters except the brightness
        parameter to generate a list of all possible combinations.
        """
        params = (
            self.all_params_except_brightness()
            if except_brightness
            else self.all_params()
        )
        return list(itertools.product(*(param.values for param in params)))

    def has_time_param(self) -> bool:
        """
        Checks if a time parameter exists in the collection.
        """
        return self.time_param is not None

    def get_time_param(self) -> Param:
        """
        Retrieves the first time parameter from the collection.
        """
        if self.has_time_param():
            return self.time_param
        raise RuntimeError(f"Time parameter missing from parameters")

    def has_brightness_param(self) -> bool:
        """
        Checks if a brightness parameter exists in the collection.
        """
        return self.brightness_param is not None

    def get_brightness_param(self) -> Param:
        """
        Retrieves the first brightness parameter from the collection.
        """
        if self.has_brightness_param():
            return self.brightness_param
        raise RuntimeError(f"Brightness parameter missing from parameters")

    def all_names_except_brightness(self) -> List[str]:
        """
        Returns the names of all parameters except the brightness parameter.
        """
        res = list(self.other.keys())
        if self.has_time_param():
            res.append(self.time_param.name)
        return res

    def all_names_except_time(self) -> List[str]:
        """
        Returns the names of all parameters except the time parameter.
        """
        res = list(self.other.keys())
        if self.has_brightness_param():
            res.append(self.brightness_param.name)
        return res

    def all_params_except_brightness(self) -> List[Param]:
        """
        Returns all parameters except the brightness parameter.
        """
        res = list(self.other.values())
        if self.has_time_param():
            res.append(self.time_param)
        return res

    def all_params_except_time(self) -> List[Param]:
        """
        Returns all parameters except the time parameter.
        """
        res = list(self.other.values())
        if self.has_brightness_param():
            res.append(self.brightness_param)
        return res

    def other_params(self):
        return list(self.other.values())

    def other_names(self) -> List[str]:
        return list(self.other.keys())

    def all_params(self) -> List[Param]:
        res = list(self.other.values())
        if self.has_time_param():
            res.append(self.time_param)
        if self.has_brightness_param():
            res.append(self.brightness_param)
        return res

    def all_names(self) -> List[str]:
        res = list(self.other.keys())
        if self.has_time_param():
            res.append(self.time_param.name)
        if self.has_brightness_param():
            res.append(self.brightness_param.name)
        return res

    def merge(self, other: "Params"):
        """
        Merges another Params object into this one.
        Raises an error if both contain a time or brightness parameter.
        """
        # Check for time param conflict
        if self.time_param and other.time_param:
            raise RuntimeError(
                f"Cannot merge: both Params instances define a time parameter "
                f"('{self.time_param.name}' and '{other.time_param.name}')"
            )
        if other.time_param:
            self.time_param = other.time_param

        # Check for brightness param conflict
        if self.brightness_param and other.brightness_param:
            raise RuntimeError(
                f"Cannot merge: both Params instances define a brightness parameter "
                f"('{self.brightness_param.name}' and '{other.brightness_param.name}')"
            )
        if other.brightness_param:
            self.brightness_param = other.brightness_param

        # Merge "other" parameters
        for param in other.other.values():
            self.add(param)

    def __str__(self):
        all_params = self.all_params()
        out = f"Params list (length {len(all_params)}):"
        for param in all_params:
            out += f"\n• {param}"
        return out

    def __eq__(self, other):
        if not isinstance(other, Params):
            self.logger.warning("Comparison failed: other is not a Params instance")
            return False

        mismatched = []

        def check(p1: Optional[Param], p2: Optional[Param], name):
            if p1 != p2:
                mismatched.append(name)
                return False
            return True

        equal = True

        equal &= check(self.time_param, other.time_param, "time_param")
        equal &= check(
            self.brightness_param, other.brightness_param, "brightness_param"
        )

        all_keys = set(self.other.keys()).union(other.other.keys())
        for key in all_keys:
            p1 = self.other.get(key)
            p2 = other.other.get(key)
            equal &= check(p1, p2, key)

        if mismatched:
            self.logger.body(
                f'Params mismatch in: {", ".join(mismatched)}',
            )
        return equal


def parse_config_param(
    name: str,
    info: Dict,
    param_type: ParamType = ParamType.OTHER,
) -> Param:
    """
    Generate a list of possible values for the parameter using settings from the config file.

    :param_name: Name of parameter as in config file.
    :param_info: Dictionary corresponding to the JSON data under the given parameter in the config file.
    :param param_type: ParamType indicating whether the Param is related to time, brightness, or neither.
    """
    logger = CustomLogger()
    logger.subheader(f"Parsing config parameter '{name}'", newline=True)

    if info["type"] == "list":
        res = ListParam(
            name,
            info["list"],
            param_type=param_type,
        )

    elif info["type"] == "range":
        res = RangeParam(
            name,
            info["range"]["minval"],
            info["range"]["maxval"],
            info["range"]["step"],
            param_type=param_type,
        )

    elif info["type"] == "logrange":
        res = LogRangeParam(
            name,
            info["logrange"]["minval"],
            info["logrange"]["maxval"],
            info["logrange"]["base"],
            info["logrange"]["n"],
            n_digits=info["logrange"]["n_digits"],
            param_type=param_type,
        )

    elif info["type"] == "random":
        res = RandomParam(
            name,
            info["random"]["minval"],
            info["random"]["maxval"],
            info["random"]["n"],
            n_digits=info["random"]["n_digits"],
            param_type=param_type,
        )

    elif info["type"] == "random_inrange":
        res = RandomInRangeParam(
            name,
            info["random_inrange"]["valid_ranges"],
            info["random_inrange"]["n"],
            param_type=param_type,
        )

    else:
        raise RuntimeError(
            "Type must be one of the following: list, range, logrange, random, random_inrange."
        )

    logger.body(f"Result: {res.__str__()}")
    return res


def parse_config_params(
    model_settings: Dict,
    time_param_name: str = "peak_mjd",
    brightness_param_name: str = "peak_appmag",
) -> Params:
    """
    Parse the parameters in the config file and generate lists of possible values for each parameter.
    """
    logger = CustomLogger()

    params: Params = Params()
    for param_name in model_settings["parameters"]:
        if param_name == time_param_name:
            param_type = ParamType.TIME
        elif param_name == brightness_param_name:
            param_type = ParamType.BRIGHTNESS
        else:
            param_type = ParamType.OTHER

        param = parse_config_param(
            param_name,
            model_settings["parameters"][param_name],
            param_type=param_type,
        )
        params.add(param)
    params.validate()
    logger.body(f"{params}", newline=True)
    return params


class SimTable(pdastrostatsclass):
    def __init__(self, brightness: float, **kwargs):
        """
        Initialize a SimTable.

        :brightness: Brightness (e.g., peak apparent magnitude or flux) for all simulations in this table.
        """
        pdastrostatsclass.__init__(self, **kwargs)
        self.logger = CustomLogger(self.__class__.__name__)
        self.brightness = brightness

    def add_row(self, data: Dict):
        """
        Add a row to the end of the table.

        :param data: Dictionary of column-value pairs.
        """
        self.newrow(data)
        # self.t = pd.concat([self.t, pd.DataFrame([data])], ignore_index=True)

    def get_sim_filename(self, model_name, tables_dir):
        return (
            f"{tables_dir}/sim_{model_name}_{format_float_string(self.brightness)}.txt"
        )

    def save_sim_table(self, model_name, tables_dir, verbose=False):
        filename = self.get_sim_filename(model_name, tables_dir)
        if verbose:
            self.logger.saving(f"Saving SimTable {filename}")
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
        self.logger = CustomLogger(self.__class__.__name__)
        self.d: Dict[str, SimTable] = {}
        self.brightness_param: Optional[Param] = None
        self.model_name = model_name

    def set_brightness_param(self, brightness_param: Param):
        if not brightness_param.is_brightness_param:
            raise ValueError(
                f"Brightness parameter must be of ParamType.BRIGHTNESS (got {brightness_param.param_type})"
            )
        if not brightness_param.has_values():
            raise ValueError(f"Brightness parameter must contain values")

        self.brightness_param = brightness_param

    def generate(
        self,
        params: Params,
        filename: Optional[str] = None,
        mjd_colname: Optional[bool] = False,
        mag_colname: Optional[bool] = False,
        flux_colname: Optional[bool] = False,
    ):
        """
        Generate the table content of each SimTable using parsed parameters from the config file.

        :param params: A Params object containing a dictionary of parameter names and their corresponding Param objects.
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
        num_rows = params.get_num_combinations()
        self.set_brightness_param(params.get_brightness_param())
        if self.brightness_param is None:
            raise RuntimeError(
                "Cannot generate SimTables: missing brightness parameter"
            )

        for brightness in self.brightness_param.values:
            self.logger.step(
                f"Generating {num_rows}-length SimTable for {self.brightness_param.name}={brightness}",
                newline=False,
            )
            self.d[brightness] = SimTable(brightness)

            combinations = params.get_combinations()
            combinations_dicts = [
                dict(zip(params.all_names_except_brightness(), combo))
                for combo in combinations
            ]

            for combination in combinations_dicts:
                combination.update({self.brightness_param.name: brightness})
                combination.update(row)
                self.d[brightness].add_row(combination)

        self.logger.success()

    def save_all(self, tables_dir: str):
        self.logger.saving(f"Saving SimTables in directory: {tables_dir}", newline=True)

        if self.brightness_param is None:
            raise RuntimeError("Cannot save SimTables: missing brightness parameter")

        if not self.d:
            raise RuntimeError(
                "Cannot save SimTables: empty dictionary (call generate() first)"
            )

        make_dir_if_not_exists(tables_dir)

        for brightness in self.brightness_param.values:
            self.d[brightness].save_sim_table(self.model_name, tables_dir)

    def load_all(self, tables_dir: str, brightness_param: Param):
        self.logger.loading(
            f"Loading SimTables in directory: {tables_dir}", newline=True
        )
        self.d = {}
        self.set_brightness_param(brightness_param)
        if self.brightness_param is None:
            raise RuntimeError("Cannot load SimTables: missing brightness parameter")

        for brightness in self.brightness_param.values:
            self.d[brightness] = SimTable(brightness)
            self.d[brightness].load_sim_table(self.model_name, tables_dir)
        self.logger.success()


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


def get_model_settings(model_name: str, step1_config: Dict) -> Dict:
    if " " in model_name:
        raise RuntimeError("Model name cannot have spaces.")
    if model_name not in step1_config:
        raise RuntimeError(
            f"Model '{model_name}' not found in simulation config file\n"
            f"Available models: {', '.join(step1_config.keys())}"
        )
    CustomLogger.s_loading(f"Getting settings for model '{model_name}'", newline=True)
    res = step1_config[model_name]
    CustomLogger.s_success()
    return res


if __name__ == "__main__":
    logger = CustomLogger()

    args = define_args().parse_args()
    config = load_config(args.config_file)
    step1_config = load_json_config(args.step1_config_file)
    model_settings = get_model_settings(args.model_name, step1_config)

    logger.header("Parsing model parameters")
    params = parse_config_params(
        model_settings,
        time_param_name=model_settings["time_parameter_name"],
        brightness_param_name=model_settings["brightness_parameter_name"],
    )
    filename, mjd_colname, mag_colname, flux_colname = parse_colname_info(
        model_settings, args.model_name
    )

    logger.header("Generating SimTables")
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
