#!/usr/bin/env python

from abc import ABC, abstractmethod
import bisect
from configparser import ConfigParser
import configparser
from functools import reduce
from getpass import getpass
from typing import Callable, Dict, Any, List, Optional, Self, Set, Tuple, Type
import re, json, requests, time, sys, io, os
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time
from collections import OrderedDict
from pdastro import pdastrostatsclass
import numpy as np
import pandas as pd
from copy import deepcopy
from pathlib import Path


# number of days to subtract from TNS discovery date to make sure no SN flux before discovery date
DISC_DATE_BUFFER = 20

# ATLAS template change dates
TEMPLATE_CHANGE_1_MJD = 58417
TEMPLATE_CHANGE_2_MJD = 58882

CONFIG_CUT_NAMES = ["uncert_cut", "x2_cut", "controls_cut", "badday_cut", "averaging"]

ATLAS_API_COLUMN_NAMES = [
    "MJD",
    "m",
    "dm",
    "uJy",
    "duJy",
    "F",
    "err",
    "chi/N",
    "RA",
    "Dec",
    "x",
    "y",
    "maj",
    "min",
    "phi",
    "apfit",
    "Sky",
    "ZP",
    "Obs",
    "Mask",
]


def add_static_methods(cls):
    """
    Class decorator that automatically adds static versions of all public instance methods.

    For each instance method (not starting with "_" or "s_"), this decorator adds a static method
    with the same name prefixed by 's_'. The static version will instantiate the class with default
    parameters and call the corresponding instance method.

    This is best used when instance methods don't depend on constructor arguments.

    Example:
        @add_static_methods
        class CustomLogger:
            def step(self, message):
                print(f"⚙️ {message}")

        CustomLogger.s_step("This works statically!")
    """
    default_instance = cls()

    for name, method in list(cls.__dict__.items()):
        if callable(method) and not name.startswith("_") and not name.startswith("s_"):

            def make_static(meth_name):
                def static_method(*args, **kwargs):
                    return getattr(default_instance, meth_name)(*args, **kwargs)

                return static_method

            setattr(cls, f"s_{name}", staticmethod(make_static(name)))

    return cls


@add_static_methods
class CustomLogger:
    def __init__(self, prefix=""):
        self.prefix = f"[{prefix}] " if prefix else ""

    def _print(
        self, message: str, symbol: str = "", newline: bool = False, dots: bool = False
    ):
        newline_part = "\n" if newline else ""
        suffix = "..." if dots else ""

        # capitalize first letter of message
        if message:
            message_part = message[0].upper() + message[1:]

        symbol_part = f"{symbol} " if symbol else ""

        print(f"{newline_part}{symbol_part}{self.prefix}{message_part}{suffix}")

    def warning(self, message: str, newline: bool = False, dots: bool = False):
        self._print(message, symbol="⚠️  WARNING:", newline=newline, dots=dots)

    def error(self, message: str, newline: bool = False, dots: bool = False):
        self._print(message, symbol="❌ ERROR:", newline=newline, dots=dots)

    def success(
        self, message: str = "Success", newline: bool = False, dots: bool = False
    ):
        self._print(message, symbol="✅", newline=newline, dots=dots)

    def header(
        self,
        message: str,
        num_dashes: int = 3,
        newline: bool = True,
    ):
        prefix_newline = "\n" if newline else ""
        dashes = "-" * num_dashes
        # no prefix for header
        print(f"{prefix_newline}{dashes} {message} {dashes}")

    def subheader(self, message: str, newline: bool = True):
        self.header(message, num_dashes=2, newline=newline)

    def step(self, message: str, newline: bool = True, dots: bool = False):
        self._print(message, symbol="⚙️ ", newline=newline, dots=dots)

    def body(self, message: str, newline: bool = False, dots: bool = False):
        self._print(message, newline=newline, dots=dots)

    def listitem(self, message: str, symbol="•", dots: bool = False):
        # no prefix for listitem
        print(f'{symbol} {message}{"..." if dots else ""}')

    def info(self, message: str, newline: bool = False, dots: bool = False):
        self._print(message, symbol="🔧", newline=newline, dots=dots)

    def secret(self, message: str, newline: bool = False, dots: bool = False):
        self._print(message, symbol="🔒", newline=newline, dots=dots)

    def loading(self, message: str, newline: bool = False, dots: bool = False):
        self._print(message, symbol="🔄", newline=newline, dots=dots)

    def saving(self, message: str, newline: bool = False, dots: bool = False):
        self._print(message, symbol="💾", newline=newline, dots=dots)

    def api(self, message: str, newline: bool = False, dots: bool = False):
        self._print(message, symbol="🌐", newline=newline, dots=dots)

    def plot(self, message: str, newline: bool = False, dots: bool = False):
        self._print(message, symbol="📈", newline=newline, dots=dots)


# convert flux to magnitude
def flux2mag(flux: float):
    """
    Convert flux in microjanskeys to apparent magnitude
    """
    return -2.5 * np.log10(flux) + 23.9


# convert magnitude to flux
def mag2flux(mag: float):
    """
    Convert apparent magnitude to flux in microjanskeys.
    """
    return 10 ** ((mag - 23.9) / -2.5)


def mag2count(mag: float, zpt: float = 20.44):
    """
    Convert apparent magnitude to TESS counts per second.
    """
    return 10 ** ((zpt - mag) / 2.5)


def count2mag(count: float, zpt: float = 20.44):
    """
    Convert TESS counts per second to apparent magnitude.
    """
    return -2.5 * np.log10(count) + zpt


def AandB(A, B) -> List:
    return list(np.intersect1d(A, B, assume_unique=False))


def AnotB(A, B) -> List:
    return list(np.setdiff1d(A, B))


def AorB(A, B) -> List:
    return list(np.union1d(A, B))


def not_AandB(A, B) -> List:
    return list(np.setxor1d(A, B))


# print iterations progress
# from https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
def print_progress_bar(
    iteration,
    total,
    prefix="",
    suffix="",
    decimals=1,
    length=100,
    fill="█",
    printEnd="\r",
):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + "-" * (length - filledLength)
    print(f"\r{prefix} |{bar}| {percent}% {suffix}", end=printEnd)
    if iteration == total:
        print()


def app2absmag(values: List[float], distance_modulus=29.04, precision=2):
    """
    Convert a list of apparent magnitude values to absolute magnitude values.

    Parameters:
    - values: list or array of apparent magnitude values
    - distance_modulus: float, the distance modulus (default 29.04)
    - precision: int, number of decimal places to round to

    Returns:
    - list of converted absolute magnitude values (as floats)
    """
    return [round(v - distance_modulus, precision) for v in values]


def load_json_config(filename: str):
    try:
        CustomLogger.s_loading(f"Loading JSON config file at {filename}", newline=True)
        with open(filename) as cfg:
            res = json.load(cfg)
            CustomLogger.s_success()
            return res
    except Exception as e:
        raise RuntimeError(f"Could not load JSON config file at {filename}: {str(e)}")


def nan_if_none(x):
    return x if x is not None else np.nan


def new_row(t: Optional[pd.DataFrame], d: Optional[Dict] = None):
    if d is None:
        d = {}

    new_row_df = pd.DataFrame([d])
    if t is None or t.empty:
        t = new_row_df
    else:
        new_row_df = new_row_df.reindex(columns=t.columns)
        t = pd.concat([t, new_row_df], axis=0, ignore_index=True)
    return t


def abbreviate_list(l: List, max_length: int = 30, abbrev_length: int = 10) -> str:
    if len(l) > max_length:
        half_abbrev_length = int(abbrev_length / 2)
        return (
            "["
            + ", ".join(map(str, l[:half_abbrev_length]))
            + ", ..., "
            + ", ".join(map(str, l[-half_abbrev_length:]))
            + f"] (length: {len(l)})"
        )
    else:
        return str(l) + f" (length: {len(l)})"


def get_allowed_presets(config: ConfigParser) -> list[str]:
    """
    Extract all preset names from the config that match 'column_name_preset.<PRESET NAME>'.
    """
    pattern = re.compile(r"^column_name_preset\.(.+)$")
    return [match.group(1) for key in config.keys() if (match := pattern.match(key))]


def parse_config_str(value: str | None):
    """Parse value from config file by converting string to string, None, or boolean."""
    if value:
        stripped = value.strip().lower()
        if stripped == "none":
            return None
        if stripped == "true":
            return True
        if stripped == "false":
            return False
    return value


def parse_comma_separated_string(string: Optional[str]):
    if string is None:
        return None

    try:
        return [item.strip() for item in string.split(",")]
    except Exception as e:
        raise RuntimeError(
            f"Could not parse comma-separated string: {string}" f"\nERROR: {str(e)}"
        )


def make_dir_if_not_exists(directory):
    """
    Creates a directory if it does not exist. Handles permission errors and other exceptions.

    :param directory: Path to the directory to create.
    """
    if not os.path.isdir(directory):
        try:
            os.makedirs(directory)
        except PermissionError:
            raise PermissionError(
                f"Permission denied: Cannot create directory at {directory}"
            )
        except FileExistsError:
            # This can occur if the directory is created between the `isdir` check and `makedirs` call.
            raise FileExistsError(f"Directory already exists: {directory}")
        except Exception as e:
            raise RuntimeError(
                f"An error occurred while creating directory {directory}: {str(e)}"
            )


# load a .ini config file
def load_config(filename):
    cfg = configparser.ConfigParser()
    try:
        CustomLogger.s_loading(f"Loading config file at {filename}", newline=True)
        cfg.read(filename)
    except Exception as e:
        raise RuntimeError(f"Could not load config file at {filename}: {str(e)}")
    CustomLogger.s_success()
    return cfg


def extract_from_subdir(
    directory: str,
    pattern: re.Pattern,
    group_name: str | int,
    convert_function: Callable = lambda x: x,
):
    if not os.path.isdir(directory):
        CustomLogger.s_warning(
            f"Cannot search because the path does not exist: {directory}"
        )
        return []

    extracted_values = set()

    for file in os.listdir(directory):
        match = pattern.match(file)
        if match:
            value = match.group(group_name)
            extracted_values.add(convert_function(value))

    if not extracted_values:
        CustomLogger.s_warning(
            f"Could not find {group_name} from the files in {directory}"
        )

    return list(extracted_values)


def has_match(directory: str, pattern: re.Pattern):
    if not os.path.isdir(directory):
        return False

    for file in os.listdir(directory):
        if pattern.match(file):
            return True
    return False


def find_all_filts(directory: str, tnsname: str) -> List[str]:
    pattern = re.compile(
        rf"^{re.escape(tnsname)}"  # tnsname
        r"(?:_i\d{3})?"  # optional control index
        r"\.(?P<filt>\w+)"  # filter (captured)
        r"(?:\.\d+\.\d+days)?"  # optional mjdbinsize
        r"(?:\.clean)?"  # optional 'clean'
        r"\.lc\.txt$"  # ends with '.lc.txt'
    )
    subdir = os.path.join(directory, tnsname)
    return extract_from_subdir(subdir, pattern, "filt")


def check_filts_against_preset(
    preset: str, allowed_presets: List[str], filts: List[str] | str
):
    if isinstance(filts, str):
        filts = [filts]

    for filt in filts:
        if filt in allowed_presets and filt != preset:
            ans = (
                input(
                    f"WARNING: Filter '{filt}' is already defined in the config file as a preset, but does not match the current preset '{preset}'. Consider removing the filter from the current list of filters to clean, then running a separate command using the preset '{filt}'. \nCONTINUE (not recommended)? (y/n): "
                )
                .strip()
                .lower()
            )
            if ans not in ["y", "yes"]:
                sys.exit(0)


def find_all_control_indices(directory: str, tnsname: str, filt=None) -> List:
    if filt is None:
        filt_pattern = r".*"
    else:
        filt_pattern = re.escape(filt)

    pattern = re.compile(
        rf"^{re.escape(tnsname)}_i(?P<control_index>\d{{3}})"  # captures control index
        rf"\.{filt_pattern}"  # match specific filt if provided
        r"(?:\.\d+\.\d+days)?"  # optional mjdbinsize
        r"(?:\.clean)?"  # optional 'clean'
        r"\.lc\.txt$"  # ends with '.lc.txt'
    )
    subdir = os.path.join(directory, tnsname, "controls")
    return extract_from_subdir(subdir, pattern, "control_index", convert_function=int)


def is_sn_in_subdir(directory: str, tnsname: str) -> bool:
    pattern = re.compile(
        rf"^{re.escape(tnsname)}"  # tnsname
        r"\.[^\.]+"  # filter
        r"(?:\.clean)?"  # optional '.clean'
        r"(?:\.\d+\.\d+days)?"  # optional mjdbinsize
        r"\.lc\.txt$"  # ends with '.lc.txt'
    )
    subdir = os.path.join(directory, tnsname)
    return has_match(subdir, pattern)


def validate_mjd_ranges(
    ranges: List[List[float]], var_name: str = "MJD_RANGES"
) -> None:
    """
    Validates that a list of MJD ranges is properly formatted.
    """
    if not isinstance(ranges, list):
        raise TypeError(f"{var_name} must be a list, got {type(ranges).__name__}")

    for i, r in enumerate(ranges):
        if not (isinstance(r, list) and len(r) == 2):
            raise TypeError(f"{var_name}[{i}] must be a 2-element list, got: {r}")
        if not all(isinstance(x, (int, float, np.integer, np.floating)) for x in r):
            raise TypeError(f"{var_name}[{i}] must contain only integers, got: {r}")
        if r[0] > r[1]:
            raise ValueError(f"{var_name}[{i}] has start > end: {r}")


def _merge_ranges(ranges: List[List[float]]) -> List[List[float]]:
    """
    Merges a list of [start, end] MJD ranges that may overlap or be adjacent.
    """
    if not ranges:
        return []

    ranges.sort()
    merged = [ranges[0]]

    for start, end in ranges[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1][1] = max(last_end, end)
        else:
            merged.append([start, end])

    return merged


def _expand_ranges(
    ranges: List[List[float]],
    min_mjd: int,
    max_mjd: int,
    expand_edges: Optional[float] = 0.0,
):
    if expand_edges == 0.0 or expand_edges is None:
        return ranges

    CustomLogger.s_body(f"Expanding range edges by {expand_edges}")

    expanded = []
    for start, end in ranges:
        new_start = max(min_mjd, start - expand_edges)
        new_end = min(max_mjd, end + expand_edges)
        expanded.append([new_start, new_end])
    if not expanded:
        return []

    # merge overlapping or adjacent ranges
    return _merge_ranges(expanded)


def get_inverse_mjd_ranges(
    mjd_ranges: List[List[float]],
    min_mjd: int,
    max_mjd: int,
    expand_edges: float = 0.0,
    exclude_mjd_ranges: Optional[List[List[float]]] = None,
) -> List[List[float]]:
    if expand_edges < 0:
        raise ValueError(
            f"Cannot expand edges of the inverse ranges by a negative amount {expand_edges}"
        )
    if len(mjd_ranges) < 1:
        return [[min_mjd, max_mjd]]
        # raise RuntimeError("MJD ranges must contain at least one range")

    validate_mjd_ranges(mjd_ranges)

    mjd_ranges = sorted(mjd_ranges)
    mjd_ranges[0][0] = max(min_mjd, mjd_ranges[0][0])
    mjd_ranges[-1][-1] = min(max_mjd, mjd_ranges[-1][-1])

    inverse = []
    cur = min_mjd
    for start, end in mjd_ranges:
        if start > cur:
            inverse.append([cur, start])
        cur = max(cur, end)

    # check if there's a gap at the end
    if cur < max_mjd:
        inverse.append([cur, max_mjd])

    # merge exclude_mjd_ranges with the inverse list
    if exclude_mjd_ranges is not None:
        validate_mjd_ranges(exclude_mjd_ranges, var_name="EXCLUDE_MJD_RANGES")
        CustomLogger.s_body(f"Excluding additional MJD ranges {exclude_mjd_ranges}")
        if len(exclude_mjd_ranges) > 0:
            combined = inverse + exclude_mjd_ranges
            inverse = _merge_ranges(combined)

    return _expand_ranges(inverse, min_mjd, max_mjd, expand_edges=expand_edges)


class StatParams:
    def __init__(self, statparams: Dict[str, int | float | None]):
        statparams = deepcopy(statparams)
        self.mean: float = nan_if_none(statparams["mean"])
        self.mean_err: float = nan_if_none(statparams["mean_err"])
        self.stdev: float = nan_if_none(statparams["stdev"])
        self.X2norm: float = nan_if_none(statparams["X2norm"])
        self.Nclip: int | float = nan_if_none(statparams["Nclip"])
        self.Ngood: int | float = nan_if_none(statparams["Ngood"])
        # self.Nexcluded: int | float = nan_if_none(statparams["Nexcluded"])
        self.ix_good: List[int] = list(statparams["ix_good"])
        self.ix_clip: List[int] = list(statparams["ix_clip"])

    def __str__(self):
        parts = []
        for key in [
            "mean",
            "mean_err",
            "stdev",
            "X2norm",
            "Nclip",
            "Ngood",
            "ix_good",
            "ix_clip",
        ]:
            val = getattr(self, key, None)
            if isinstance(val, float):
                parts.append(f"{key}={val:.17g}")  # Full float precision
            else:
                parts.append(f"{key}={val}")
        return f"StatParams({', '.join(parts)})"


class PlotLimits:
    def __init__(self, xlower=None, xupper=None, ylower=None, yupper=None):
        self.xlower = xlower
        self.xupper = xupper
        self.ylower = ylower
        self.yupper = yupper

    def set_lims(
        self,
        xlims: Optional[tuple[float | None, float | None]] = None,
        ylims: Optional[tuple[float | None, float | None]] = None,
    ):
        if xlims is not None:
            self.set_xlims(xlims)
        if ylims is not None:
            self.set_ylims(ylims)

    def set_xlims(self, xlims: tuple[float | None, float | None] | None):
        if xlims is None:
            return

        if len(xlims) != 2:
            raise ValueError(f"xlims must be a tuple of length 2, got {len(xlims)}")

        if xlims[0] is not None and xlims[1] is not None and xlims[0] >= xlims[1]:
            raise ValueError(
                f"xlims lower limit {xlims[0]} must be less than upper limit {xlims[1]}"
            )

        if xlims[0] is not None:
            self.xlower = xlims[0]
        if xlims[1] is not None:
            self.xupper = xlims[1]

    def set_ylims(self, ylims: tuple[float | None, float | None] | None):
        if ylims is None:
            return

        if len(ylims) != 2:
            raise ValueError(f"ylims must be a tuple of length 2, got {len(ylims)}")

        if ylims[0] is not None and ylims[1] is not None and ylims[0] >= ylims[1]:
            raise ValueError(
                f"ylims lower limit {ylims[0]} must be less than upper limit {ylims[1]}"
            )

        if ylims[0] is not None:
            self.ylower = ylims[0]
        if ylims[1] is not None:
            self.yupper = ylims[1]

    def get_xlims(self) -> tuple[float | None, float | None] | None:
        if self.xlower is None and self.xupper is None:
            return None
        return self.xlower, self.xupper

    def get_ylims(self) -> tuple[float | None, float | None] | None:
        if self.ylower is None and self.yupper is None:
            return None
        return self.ylower, self.yupper

    def is_empty(self):
        return (
            self.xlower is None
            and self.xupper is None
            and self.ylower is None
            and self.yupper is None
        )

    def __str__(self):
        return f"Plot limits: x-axis [{self.xlower}, {self.xupper}], y-axis [{self.ylower}, {self.yupper}]"


class FomLimits:
    def __init__(self):
        self._values = {}
        self.logger = CustomLogger(self.__class__.__name__)

    def add(self, sigma_kern: float, values: List[float | int] | float | int):
        """
        Add new FOM limit(s) to a sigma_kern.
        """
        if isinstance(values, list):
            sorted_list_to_add = self._validate_flat_list(values)
        else:
            sorted_list_to_add = [float(values)]

        if self.has(sigma_kern):
            existing = set(self._values[sigma_kern])
            for v in sorted_list_to_add:
                if v not in existing:
                    bisect.insort(self._values[sigma_kern], v)
        else:
            self._values[sigma_kern] = sorted_list_to_add

    def set(
        self,
        sigma_kerns: List[float],
        values: (
            List[float]
            | List[List[float]]
            | Dict[float, float]
            | Dict[float, List[float]]
        ),
    ):
        """
        Overwrite any current sigma_kerns and FOM limits with new ones.
        """
        if self._values:
            self.logger.warning("Overwriting current FOM limits with new ones")
        self._values = self.validate(sigma_kerns, values)

    def set_blank(self, sigma_kerns: List[float]):
        """
        Overwrite any current sigma_kerns with new ones, and any current FOM limits with 0.0.
        """
        self.set(sorted(sigma_kerns), [0.0] * len(sigma_kerns))

    def has(self, sigma_kern) -> bool:
        return sigma_kern in self._values.keys()

    def get(self, sigma_kern: float, index: int = -1) -> float:
        """
        Get an FOM limit at a certain index for a specific sigma_kern.
        """
        if not self.has(sigma_kern):
            raise ValueError(
                f"FOM limits for sigma_kern {sigma_kern} not found (existing sigma_kerns: {self._values.keys()})"
            )
        if index >= len(self._values[sigma_kern]):
            raise ValueError(
                f"Cannot get FOM limit at index {index}; only {len(self._values[sigma_kern])} FOM limits exist per sigma_kern"
            )
        return self._values[sigma_kern][index]

    def get_multi(self, sigma_kern: float) -> List[float]:
        """
        Get all FOM limits for a specific sigma_kern.
        """
        if not self.has(sigma_kern):
            raise ValueError(
                f"FOM limits for sigma_kern {sigma_kern} not found (existing sigma_kerns: {self._values.keys()})"
            )
        return self._values[sigma_kern]

    def get_all_multi(self) -> Dict[float, List[float]]:
        """
        Get all FOM limits for each sigma_kern.
        """
        return self._values

    def get_all_single(self, index: int = -1) -> Dict[float, float]:
        """
        Get an FOM limit at a certain index for each sigma_kern.
        """
        return {
            sigma_kern: self.get(sigma_kern, index=index)
            for sigma_kern in self._values.keys()
        }

    def _validate_flat_list(self, data: List | float) -> List[float]:
        """
        Verifies that all entries in a list are numbers and converts them to floats. Sorts the resulting list.

        :param data: List of numbers (int or float)
        :return: Sorted list of floats
        :raises TypeError: If any item is not a number
        """
        if isinstance(data, list):
            for i, item in enumerate(data):
                if not isinstance(item, (int, float)):
                    raise TypeError(f"Item at index {i} is not a number: {item}")
            res = [float(x) for x in data]
            res.sort()
            return res
        else:
            return [float(data)]

    def validate(
        self,
        sigma_kerns: List[float],
        fom_limits: (
            List[float]
            | List[List[float]]
            | Dict[float, float]
            | Dict[float, List[float]]
        ),
    ) -> Dict[float, List[float]]:
        """
        Validate and convert FOM limits into a dictionary with sigma_kerns as keys.
        WARNING: If passing fom_limits as a list, both fom_limits and sigma_kerns must be sorted.

        :param fom_limits: FOM limits as a list or dictionary.

            Supports:

            - List[float]: one FOM limit per sigma_kern
            - List[List[float]]: multiple FOM limits per sigma_kern
            - Dict[float, float]: one FOM limit per sigma_kern
            - Dict[float, List[float]]: multiple FOM limits per sigma_kern, already structured

        :param sigma_kerns: Kernel sizes corresponding to the FOM limits.
        """
        # List[float] or List[List[float]]
        if isinstance(fom_limits, list):
            if not fom_limits:
                raise ValueError("FOM limits list is empty")

            if len(fom_limits) != len(sigma_kerns):
                raise ValueError(
                    f"Length mismatch: got {len(fom_limits)} FOM limits but {len(sigma_kerns)} sigma_kerns"
                )

            # List[float]
            if all(isinstance(v, (int, float)) for v in fom_limits):
                # wrap each float in a list
                return dict(
                    zip(
                        sigma_kerns, [[v] for v in self._validate_flat_list(fom_limits)]
                    )
                )

            # List[List[float]]
            elif all(isinstance(v, list) for v in fom_limits):
                return {
                    sigma_kern: self._validate_flat_list(sublist)
                    for sigma_kern, sublist in zip(sigma_kerns, fom_limits)
                }

            else:
                raise TypeError(
                    "If passing a list, it must contain only numbers or only sublists of numbers"
                )

        # Dict[float, float] or Dict[float, List[float]]
        elif isinstance(fom_limits, dict):
            if set(fom_limits.keys()) != set(sigma_kerns):
                raise ValueError("FOM limits dict keys must exactly match sigma_kerns")

            return {
                sigma_kern: self._validate_flat_list(sublist)
                for sigma_kern, sublist in fom_limits.items()
            }

        else:
            raise TypeError(
                "FOM limits must be a list of floats/sublists or a dict with float/list values"
            )

    def merge(self, other: Self):
        if not isinstance(other, FomLimits):
            raise TypeError("Argument must be an instance of FomLimits")

        for other_key, other_values in other.get_all_multi().items():
            self.add(other_key, other_values)

    def __str__(self):
        return self.get_all_multi().__str__()


class PresetColumnNames:
    """
    Class to handle loading and managing column names for light curve conversion
    to ATClean-readable format, as defined in the config file.
    """

    def __init__(self, config: ConfigParser, preset: str):
        self.preset = preset
        self._read_config(config)

    def _validate_columns_dict(self, columns_dict: Dict, no_nones: bool = False):
        for key, name in columns_dict.items():
            name = parse_config_str(name)
            if no_nones and name is None:
                raise RuntimeError(
                    f"Column name '{name}' in config preset {self.preset} cannot be None"
                )
            columns_dict[key] = name
        return columns_dict

    def _read_config(self, config: ConfigParser):
        try:
            config_preset_settings = dict(config[f"column_name_preset.{self.preset}"])
        except:
            raise RuntimeError(
                f"Preset '{self.preset}' (field '{f'column_name_preset.{self.preset}'}') not found in config file."
            )

        self.required_columns: Dict[str, str] = self._validate_columns_dict(
            {
                "mjd": config_preset_settings.get("mjd_column_name"),
                "flux": config_preset_settings.get("flux_column_name"),
                "dflux": config_preset_settings.get("dflux_column_name"),
                "mjdbin": config["column_name_preset"]["mjd_bin_column_name"],
                "mask": config["column_name_preset"]["mask_column_name"],
                "fdf": config["column_name_preset"]["snr_column_name"],
            },
            no_nones=True,
        )

        self.optional_columns: Dict[str, str | None] = self._validate_columns_dict(
            {
                "chisquare": config_preset_settings.get("chisquare_column_name"),
                "filt": config_preset_settings.get("filter_column_name"),
                "mag": config_preset_settings.get("mag_column_name"),
                "dmag": config_preset_settings.get("dmag_column_name"),
                "ra": config_preset_settings.get("ra_column_name"),
                "dec": config_preset_settings.get("dec_column_name"),
                "zpt": config_preset_settings.get("zpt_column_name"),
            }
        )

        # extra columns to copy
        extra_columns = parse_config_str(config_preset_settings["extra_columns"])
        if not isinstance(extra_columns, str):
            raise ValueError(
                f"Extra columns in config file ({config_preset_settings['extra_columns']}) must be comma-separated list (got {extra_columns})"
            )
        self.extra_columns: List[str] = (
            []
            if extra_columns is None
            else [col.strip() for col in extra_columns.split(",")]
        )

    def add(
        self, key: str, name: str, is_required: bool = False, overwrite: bool = False
    ):
        if not isinstance(key, str) or not key.strip():
            raise ValueError("Column key must be a non-empty string.")
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"Column name for key '{key}' must be a non-empty string.")

        # if key exists but the name is different, raise an error
        existing_name = self.required_columns.get(key) or self.optional_columns.get(key)
        if existing_name and existing_name != name and not overwrite:
            raise RuntimeError(
                f"Column key '{key}' is already defined with a different name '{existing_name}'."
            )

        if is_required:
            self.required_columns[key] = name
        else:
            self.optional_columns[key] = name

    def add_many(self, coldict: Dict, is_required: bool = False):
        for key, name in coldict.items():
            self.add(key, name, is_required=is_required)

    def update(self, key: str, name: str, is_required: bool = False):
        if is_required:
            if key not in self.required_columns:
                raise RuntimeError(
                    f"Cannot update non-existing required column name {key} with '{name}'"
                )
            self.required_columns[key] = name
        else:
            if key not in self.optional_columns:
                raise RuntimeError(
                    f"Cannot update non-existing optional column name {key} with '{name}'"
                )
            self.optional_columns[key] = name

    def update_many(self, coldict: Dict, is_required: bool = False):
        for key, name in coldict.items():
            self.update(key, name, is_required=is_required)

    def remove(self, key: str):
        if key in self.optional_columns:
            del self.optional_columns[key]
        if key in self.optional_columns:
            del self.optional_columns[key]

    def get_required_column_names(self, is_averaged: bool = False):
        if is_averaged:
            return [
                self.required_columns["mjdbin"],
                self.required_columns["flux"],
                self.required_columns["dflux"],
                self.required_columns["mask"],
            ]

        return [
            self.required_columns["mjd"],
            self.required_columns["flux"],
            self.required_columns["dflux"],
        ]

    def get_optional_column_names(self):
        return list(self.optional_columns.values())

    def get_all_columns_to_copy(self) -> List[str]:
        """Return all columns that should be copied into the output light curve."""
        colset: Set[str] = (
            set(self.required_columns.values())
            | set(filter(None, self.optional_columns.values()))
            | set(self.extra_columns)
        )
        return list(colset)

    def has(self, name: str):
        return name in self.required_columns or name in self.optional_columns

    def __getattr__(self, name: str) -> str | None:
        """
        Dynamic access to column names, e.g., obj.mjd or obj.chisquare.
        """
        if (
            "required_columns" not in self.__dict__
            or "optional_columns" not in self.__dict__
        ):
            raise AttributeError(
                f"'{type(self).__name__}' object has no attribute '{name}'"
            )

        if name in self.required_columns:
            return self.required_columns[name]
        if name in self.optional_columns:
            return self.optional_columns[name]
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )

    def __str__(self) -> str:
        """Readable string representation of all column names."""
        skip_colnames = ["mjdbin", "fdf", "mask"]

        lines = ["-- Required Columns --"]
        for k, v in self.required_columns.items():
            if k in skip_colnames:
                continue
            lines.append(f"{k}: {v}")

        lines.append("-- Optional Columns --")
        for k, v in self.optional_columns.items():
            lines.append(f"{k}: {v}")

        lines.append("-- Extra Columns to Copy --")
        lines.append(", ".join(self.extra_columns) if self.extra_columns else "(None)")

        return "\n".join(lines)


def load_preset_column_names_from_config(
    preset: str,
    config: ConfigParser,
    filts: Optional[List[str] | str] = None,
    verbose: bool = True,
) -> PresetColumnNames:
    """
    Load a set of preset column names from a configuration file.

    This function retrieves column name presets defined in a config file and returns
    a `PresetColumnNames` instance corresponding to the specified preset. It also performs
    validation and optionally displays status messages.

    :param preset: The name of the preset to load from the configuration file. Must be one of the allowed presets defined in the config.
    :param config: A ConfigParser object containing preset definitions (typically parsed from config.ini).
    :param filts: Filters that may be used to warn if it matches a preset but differs from `preset`. Useful for catching potential mismatches.

    Raises RuntimeError if the provided preset is None or not found among the allowed presets.

    If `filt` is specified and matches a preset name different from `preset`, a warning is printed.
    """
    if verbose:
        CustomLogger.s_loading(
            f"Loading '{preset}' preset column names from config.ini", newline=True
        )

    allowed_presets = get_allowed_presets(config)
    if preset is None or preset not in allowed_presets:
        raise RuntimeError(
            f"Please specify the preset name to load from the config file (allowed presets: {allowed_presets})"
        )

    if filts is not None:
        check_filts_against_preset(preset, allowed_presets, filts)

    colnames = PresetColumnNames(config, preset)
    if verbose:
        CustomLogger.s_success()
        print(colnames.__str__())
    return colnames


class Credentials:
    def __init__(
        self, atlas_username, atlas_password, tns_api_key, tns_id, tns_bot_name
    ):
        self.atlas_username = parse_config_str(atlas_username)
        self.atlas_password = parse_config_str(atlas_password)
        self.tns_api_key = parse_config_str(tns_api_key)
        self.tns_id = parse_config_str(tns_id)
        self.tns_bot_name = parse_config_str(tns_bot_name)

    def prompt_for_atlas_password(self):
        if self.atlas_password is None:
            self.atlas_password = getpass(prompt="Enter ATLAS password: ")

    def prompt_for_tns_creds(self):
        if self.tns_id is None:
            self.tns_id = getpass(prompt="Enter TNS ID: ")
        if self.tns_bot_name is None:
            self.tns_bot_name = getpass(prompt="Enter TNS bot name: ")
        if self.tns_api_key is None:
            self.tns_api_key = getpass(prompt="Enter TNS API key: ")


class BaseAngle(ABC):
    def __init__(self, string=None):
        self.angle = None
        if string:
            self.set_angle(string)

    def set_angle(self, string):
        """Template method calling specific parsing."""
        if self._is_nan(string):
            self.angle = None
            return

        self.angle = self._parse_angle(string)

    @staticmethod
    def _is_nan(value):
        """Check if a value is NaN, None, or an empty string."""
        return (
            value is None
            or (isinstance(value, float) and np.isnan(value))
            or (isinstance(value, str) and value.strip().lower() in ["nan", ""])
        )

    @abstractmethod
    def _parse_angle(self, string):
        """Abstract method to be implemented by subclasses to parse angles."""
        pass


class RA(BaseAngle):
    def _parse_angle(self, string):
        """Parse RA angle, using hours if ':' is present, degrees otherwise."""
        s = re.compile(":")
        if isinstance(string, str) and s.search(string):
            return Angle(string, u.hour)
        else:
            return Angle(string, u.degree)


class Dec(BaseAngle):
    def _parse_angle(self, string):
        """Parse Dec angle, always using degrees."""
        return Angle(string, u.degree)


class Coordinates:
    def __init__(self, ra: str | None = None, dec: str | None = None):
        self.ra: RA = RA(ra)
        self.dec: Dec = Dec(dec)

    def set_RA(self, ra: str):
        self.ra = RA(ra)

    def set_Dec(self, dec: str):
        self.dec = Dec(dec)

    def get_RA_str(self) -> str:
        if self._is_angle_missing(self.ra):
            return str(np.nan)
        return f"{self.ra.angle.degree:0.14f}"

    def get_Dec_str(self) -> str:
        if self._is_angle_missing(self.dec):
            return str(np.nan)
        return f"{self.dec.angle.degree:0.14f}"

    def _is_angle_missing(self, angle: BaseAngle) -> bool:
        return angle.angle is None

    def is_complete(self) -> bool:
        return not self._is_angle_missing(self.ra) and not self._is_angle_missing(
            self.dec
        )

    def is_empty(self) -> bool:
        # both RA and Dec missing
        return self._is_angle_missing(self.ra) and self._is_angle_missing(self.dec)

    def is_incomplete(self) -> bool:
        # one or both of RA and Dec missing
        return self._is_angle_missing(self.ra) or self._is_angle_missing(self.dec)

    def is_ra_present(self) -> bool:
        return not self._is_angle_missing(self.ra)

    def is_dec_present(self) -> bool:
        return not self._is_angle_missing(self.dec)

    def ra_andor_dec_present(self) -> bool:
        return not self._is_angle_missing(self.ra) or not self._is_angle_missing(
            self.dec
        )

    def get_distance(self, other: Self) -> Angle:
        if self.is_incomplete():
            raise ValueError(
                "To get distance, RA and Dec must be present in this Coordinates object"
            )
        if other.is_incomplete():
            raise ValueError(
                "To get distance, RA and Dec must be present in other Coordinates object"
            )

        c1 = SkyCoord(self.ra.angle, self.dec.angle, frame="fk5")
        c2 = SkyCoord(other.ra.angle, other.dec.angle, frame="fk5")
        return c1.separation(c2)

    def __str__(self):
        output = []
        if self.is_ra_present():
            output.append(f"RA {self.get_RA_str()}")
        if self.is_dec_present():
            output.append(f"Dec {self.get_Dec_str()}")

        if len(output) < 1:
            return f"⚠️  WARNING: Coordinates are empty and cannot be printed."
        return ", ".join(output)


# input/output table containing TNS names, RA, Dec, and MJD0
# (TODO: if MJD0=None, consider entire light curve as pre-SN light curve)
class SnInfoTable:
    def __init__(self, directory, filename=None):
        self.logger = CustomLogger(self.__class__.__name__)

        if filename is None:
            self.filename = f"{directory}/sninfo.txt"
        else:
            self.filename = f"{directory}/{filename}"

        try:
            self.logger.loading(
                f"Loading SN info table at {self.filename}", newline=True
            )
            self.t = pd.read_table(self.filename, sep="\s+")
            if not "tnsname" in self.t.columns:
                raise RuntimeError('SN info table must have a "tnsname" column.')
            self.t["ra"] = self.t["ra"].astype(str)
            self.t["dec"] = self.t["dec"].astype(str)
            if "mjd0" not in self.t.columns:
                self.t["mjd0"] = np.nan
            self.logger.success()
        except Exception:
            self.logger.body(
                f"No existing SN info table at that path; creating blank table",
                dots=True,
            )
            self.t = pd.DataFrame(
                columns=["tnsname", "ra", "dec", "mjd0", "center_ra", "center_dec"]
            )

    def get_row(self, tnsname) -> tuple[int, pd.Series] | tuple[int, None]:
        if self.t.empty:
            # raise RuntimeError(f'Error: Cannot get info for SN {tnsname}--table is empty.')
            return -1, None

        matching_ix = self.t.index[self.t["tnsname"].eq(tnsname)]
        if len(matching_ix) >= 2:
            self.logger.warning(
                f"SN info table has {len(matching_ix)} matching rows for TNS name {tnsname}. Dropping duplicate rows",
                dots=True,
            )
            self.t.drop(index=matching_ix[1:], inplace=True)
            first_ix = self.t.index[self.t["tnsname"].eq(tnsname)][0]
            return first_ix, self.t.loc[first_ix]
        elif len(matching_ix) == 1:
            return matching_ix[0], self.t.loc[matching_ix[0]]
        else:
            return -1, None

    def is_nan(self, value) -> bool:
        if value is None:
            return True
        if isinstance(value, (float, np.floating)) and np.isnan(value):
            return True
        if isinstance(value, str) and value.strip().lower() in ["nan", "", "none"]:
            return True
        return False

    def get_info(self, tnsname) -> tuple[Coordinates, Coordinates, float | None]:
        _, row = self.get_row(tnsname)
        if row is None:
            return None, None, None

        assert "ra" in row
        assert "dec" in row
        ra = None if self.is_nan(row["ra"]) else row["ra"]
        dec = None if self.is_nan(row["dec"]) else row["dec"]

        center_ra, center_dec = None, None
        if "center_ra" in row and not self.is_nan(row["center_ra"]):
            center_ra = row["center_ra"]
        if "center_dec" in row and not self.is_nan(row["center_dec"]):
            center_dec = row["center_dec"]

        if "mjd0" not in row or self.is_nan(row["mjd0"]):
            mjd0 = None
        else:
            if not isinstance(row["mjd0"], (int, float, np.integer, np.floating)):
                raise RuntimeError(f"Invalid MJD0: {row['mjd0']}")
            mjd0 = float(row["mjd0"])

        return Coordinates(ra, dec), Coordinates(center_ra, center_dec), mjd0

    def update_row_at_index(
        self,
        index,
        coords: Optional[Coordinates] = None,
        mjd0: Optional[float] = None,
        overwrite=False,
    ):
        try:
            if (
                (overwrite or np.isnan(self.t.at[index, "mjd0"]))
                and mjd0 is not None
                and not np.isnan(mjd0)
            ):
                self.t.loc[index, "mjd0"] = mjd0

            if (
                (overwrite or self.is_nan(self.t.loc[index, "ra"]))
                and not coords is None
                and coords.is_ra_present()
            ):
                self.t.loc[index, "ra"] = f"{coords.ra.angle.degree:0.14f}"

            if (
                (overwrite or self.is_nan(self.t.loc[index, "dec"]))
                and not coords is None
                and coords.is_dec_present()
            ):
                self.t.loc[index, "dec"] = f"{coords.dec.angle.degree:0.14f}"
        except Exception as e:
            raise RuntimeError(
                f"Could not update SN info table at index {index}: {str(e)}"
            )

    def add_new_row(
        self,
        tnsname,
        coords: Optional[Coordinates] = None,
        mjd0: Optional[float] = None,
    ):
        if mjd0 is None:
            mjd0 = np.nan

        ra = np.nan
        dec = np.nan
        if not coords is None:
            if coords.is_ra_present():
                ra = f"{coords.ra.angle.degree:0.14f}"
            if coords.is_dec_present():
                dec = f"{coords.dec.angle.degree:0.14f}"

        row = {"tnsname": tnsname, "ra": ra, "dec": dec, "mjd0": mjd0}
        self.t = new_row(self.t, row)

    def update_row(
        self,
        tnsname,
        coords: Optional[Coordinates] = None,
        mjd0: Optional[float] = None,
        overwrite=False,
    ):
        if self.t.empty:
            self.add_new_row(tnsname, coords, mjd0)
            return

        matching_ix = np.where(self.t["tnsname"].eq(tnsname))[0]
        if len(matching_ix) > 1:
            raise RuntimeError(
                f"SN info table has {len(matching_ix)} matching rows for TNS name {tnsname}."
            )
        elif len(matching_ix) == 1:
            index = matching_ix[0]
            self.update_row_at_index(
                index, coords=coords, mjd0=mjd0, overwrite=overwrite
            )
        else:
            self.add_new_row(tnsname, coords=coords, mjd0=mjd0)

    def save(self):
        self.logger.saving(f"Saving SN info table at {self.filename}", newline=True)
        self.t["ra"] = self.t["ra"].astype(str)
        self.t["dec"] = self.t["dec"].astype(str)
        self.t.to_string(self.filename, index=False)

    def __str__(self):
        return self.t.to_string()


def format_float_string(value: int | float):
    """
    Format with up to 5 decimals, strip trailing zeros, then ensure at least 1 decimal
    """
    formatted = f"{value:.5f}".rstrip("0").rstrip(".")
    if "." not in formatted:
        formatted += ".0"
    return formatted


def get_filepath(
    directory, tnsname, filt="o", control_index=0, mjdbinsize=None, cleaned=False
):
    filename = f"{directory}/{tnsname}"

    if control_index != 0:
        filename += "/controls"

    filename += f"/{tnsname}"

    if control_index != 0:
        filename += f"_i{control_index:03d}"

    filename += f".{filt}"

    if mjdbinsize:
        filename += f".{format_float_string(mjdbinsize)}days"

    if cleaned:
        filename += f".clean"

    filename += ".lc.txt"
    return filename


def query_tns(tnsname, api_key, tns_id, bot_name):
    logger = CustomLogger("query_tns")

    if tns_id is None or bot_name is None:
        logger.warning(
            "Cannot query TNS without TNS ID and bot name. Please specify these parameters in config.ini."
        )
        return None

    json_data = None
    try:
        url = "https://www.wis-tns.org/api/get/object"
        json_file = OrderedDict(
            [("objname", tnsname), ("objid", ""), ("photometry", "1"), ("spectra", "1")]
        )
        data = {"api_key": api_key, "data": json.dumps(json_file)}
        response = requests.post(
            url,
            data=data,
            headers={
                "User-Agent": 'tns_marker{"tns_id":"%s","type": "bot", "name":"%s"}'
                % (tns_id, bot_name)
            },
        )
        json_data = json.loads(response.text, object_pairs_hook=OrderedDict)
        return json_data
    except Exception as e:
        if json_data and "data" in json_data:
            logger.body(json_data["data"])
        else:
            logger.body("No JSON data received or failed to parse response")
        raise RuntimeError("ERROR in query_tns(): " + str(e))


def get_tns_coords_from_json(json_data):
    try:
        coords = Coordinates(json_data["data"]["ra"], json_data["data"]["dec"])
        return coords
    except Exception as e:
        raise RuntimeError(f"Failed to get coordinates from TNS JSON data: {str(e)}")


def get_tns_mjd0_from_json(json_data, use_disc_date_buffer: bool = True):
    try:
        disc_date = json_data["data"]["discoverydate"]
        date = list(disc_date.partition(" "))[0]
        time = list(disc_date.partition(" "))[2]
        date_object = Time(date + "T" + time, format="isot", scale="utc")

        mjd0 = date_object.mjd
        if use_disc_date_buffer:
            mjd0 -= DISC_DATE_BUFFER
        return mjd0
    except Exception as e:
        raise RuntimeError(f"Failed to get discovery date from TNS JSON data: {str(e)}")


def resolve_mjd0(
    tnsname: str,
    sninfo: SnInfoTable,
    credentials: Credentials,
    use_disc_date_buffer: bool = True,
) -> float | None:
    logger = CustomLogger("get_mjd0_from_tns")
    logger.subheader("Resolving MJD0")

    _, sninfo_row = sninfo.get_row(tnsname)
    if not sninfo_row is None and not np.isnan(sninfo_row["mjd0"]):
        # get MJD0 from SN info table
        logger.info(
            f'Setting MJD0 to {sninfo_row["mjd0"]} MJD from SN info table', newline=True
        )
        mjd0 = float(sninfo_row["mjd0"])
        if isinstance(mjd0, (int, float)):
            logger.success()
            return mjd0
        else:
            logger.warning(f"Cannot convert to float: {sninfo_row['mjd0']}")

    # try querying TNS
    logger.api(f"Querying TNS for SN {tnsname} discovery date", newline=True)
    credentials.prompt_for_tns_creds()
    json_data = query_tns(
        tnsname,
        credentials.tns_api_key,
        credentials.tns_id,
        credentials.tns_bot_name,
    )
    mjd0 = get_tns_mjd0_from_json(json_data, use_disc_date_buffer=use_disc_date_buffer)
    if mjd0 is not None and not np.isnan(mjd0):
        logger.success()
        return mjd0
    else:
        logger.warning(
            "Could not resolve SN MJD0 from command line, SnInfoTable, or TNS discovery date"
        )
        return None


def get_tns_data(
    tnsname: str, tns_api_key, tns_id, tns_bot_name, use_disc_date_buffer: bool = True
) -> tuple[float, Coordinates]:
    logger = CustomLogger("get_tns_data")
    logger.api(f"Querying TNS for {tnsname} RA, Dec, and MJD0", newline=True)

    json_data = query_tns(tnsname, tns_api_key, tns_id, tns_bot_name)
    if json_data is None:
        logger.warning(f"No data returned; skipping", dots=True)
        return

    coords = get_tns_coords_from_json(json_data)
    mjd0 = get_tns_mjd0_from_json(json_data, use_disc_date_buffer=use_disc_date_buffer)
    logger.success()
    return mjd0, coords


def query_atlas(headers, ra, dec, min_mjd, max_mjd):
    logger = CustomLogger("query_atlas")

    baseurl = "https://fallingstar-data.com/forcedphot"
    task_url = None
    while not task_url:
        with requests.Session() as s:
            resp = s.post(
                f"{baseurl}/queue/",
                headers=headers,
                data={
                    "ra": ra,
                    "dec": dec,
                    "send_email": False,
                    "mjd_min": min_mjd,
                    "mjd_max": max_mjd,
                },
            )
            if resp.status_code == 201:
                task_url = resp.json()["url"]
                logger.body(f"Task url: {task_url}")
            elif resp.status_code == 429:
                message = resp.json()["detail"]
                logger.body(f"{resp.status_code} {message}")
                t_sec = re.findall(r"available in (\d+) seconds", message)
                t_min = re.findall(r"available in (\d+) minutes", message)
                if t_sec:
                    waittime = int(t_sec[0])
                elif t_min:
                    waittime = int(t_min[0]) * 60
                else:
                    waittime = 10
                logger.body(f"Waiting {waittime} seconds")
                time.sleep(waittime)
            else:
                raise RuntimeError(
                    f"Error querying ATLAS API: {resp.status_code} {resp.text}"
                )
                # logger.error(f"{resp.status_code}")
                # logger.body(resp.text)
                # sys.exit()

    result_url = None
    taskstarted_printed = False

    logger.body("Waiting for job to start...")
    while not result_url:
        with requests.Session() as s:
            resp = s.get(task_url, headers=headers)
            if resp.status_code == 200:
                if not (resp.json()["finishtimestamp"] is None):
                    result_url = resp.json()["result_url"]
                    logger.success(
                        f"Task is complete with results available at {result_url}"
                    )
                    break
                elif resp.json()["starttimestamp"]:
                    if not taskstarted_printed:
                        logger.body(
                            f"Task is running (started at {resp.json()['starttimestamp']})"
                        )
                        taskstarted_printed = True
                    time.sleep(2)
                else:
                    # print(f"Waiting for job to start (queued at {resp.json()['timestamp']})")
                    time.sleep(4)
            else:
                logger.error(f"{resp.status_code}")
                logger.body(resp.text)
                sys.exit()

    with requests.Session() as s:
        if result_url is None:
            logger.warning("Empty light curve (no data within the queried MJD range)")
            dfresult = pd.DataFrame(columns=ATLAS_API_COLUMN_NAMES)
        else:
            result = s.get(result_url, headers=headers).text
            dfresult = pd.read_csv(io.StringIO(result.replace("###", "")), sep="\s+")

    return dfresult


def reformat_dir(directory: str, tnsname: str):
    """Removes tnsname and trailing slash from the end of directory path if present."""
    if directory.endswith(tnsname):
        return directory[: -len(tnsname)].rstrip(os.sep)
    return directory


def hexstring_to_int(hexstring):
    return int(hexstring, 16)


def combine_flags(flags: List[int]) -> int:
    return reduce(lambda x, y: x | y, flags, 0)


def get_config_custom_cuts(config: ConfigParser) -> List:
    logger = CustomLogger("get_config_custom_cuts")
    logger.body("Searching config file for custom cuts", newline=True)

    required_keys = {"column", "flag", "max_value", "min_value"}
    custom_cuts = []

    for key in config:
        if key.endswith("_cut") and not key in CONFIG_CUT_NAMES:
            if not required_keys.issubset(config[key].keys()):
                logger.warning(
                    f"Custom cut {key} missing required fields (required fields: {required_keys}); skipping",
                    dots=True,
                )
            else:
                custom_cuts.append(config[key])

    logger.success(f"Found {len(custom_cuts)} custom cuts")
    return custom_cuts


def get_config_flags(config: ConfigParser) -> int:
    # get main flags for Uncertainty Cut, Chi-Square Cut, Control Light Curve Cut, and Bad Day Cut
    flags = [
        hexstring_to_int(config["uncert_cut"]["flag"]),
        hexstring_to_int(config["x2_cut"]["flag"]),
        hexstring_to_int(config["controls_cut"]["bad_flag"]),
        hexstring_to_int(config["averaging"]["flag"]),
    ]

    # add custom cuts flags to flags list
    custom_cuts = get_config_custom_cuts(config)
    for cut in custom_cuts:
        flags.append(hexstring_to_int(cut["flag"]))

    return combine_flags(flags)


class Cut(ABC):
    def __init__(
        self,
        column: Optional[str] = None,
        flag: Optional[int] = None,
        min_value: Optional[float] = None,
        max_value: Optional[float] = None,
    ):
        self.column = column
        self.flag = flag

        if min_value is not None and max_value is not None and min_value >= max_value:
            raise ValueError(
                f"Min value {min_value} cannot be greater than or equal to max value {max_value}"
            )
        self.min_value = min_value
        self.max_value = max_value

    def can_apply_directly(self) -> bool:
        return (
            self.flag is not None
            and self.column is not None
            and self.column != ""
            and (self.min_value is not None or self.max_value is not None)
        )

    @staticmethod
    @abstractmethod
    def name() -> str:
        pass

    def get_flags(self, return_keys=False) -> Dict[str, int] | List[int]:
        """Extract all attributes that end in '_flag' or are 'flag'."""
        flags = {
            key: value
            for key, value in vars(self).items()
            if isinstance(value, int) and (key.endswith("_flag") or key == "flag")
        }
        return flags if return_keys else list(flags.values())

    def __str__(self) -> str:
        details = []
        if self.column:
            details.append(f"column={self.column}")
        if self.min_value is not None:
            details.append(f"min_value={self.min_value}")
        if self.max_value is not None:
            details.append(f"max_value={self.max_value}")

        flags = self.get_flags(return_keys=True)
        if flags:
            for flag_name, flag_value in flags.items():
                details.append(f"{flag_name}={hex(flag_value)}")

        return (f"{self.name()}: " + ", ".join(details)) if details else self.name()


class CustomCut(Cut):
    def __init__(
        self,
        column: str,
        flag: int,
        min_value: Optional[float] = None,
        max_value: Optional[float] = None,
    ):
        super().__init__(
            column=column, flag=flag, min_value=min_value, max_value=max_value
        )
        if not self.can_apply_directly():
            raise ValueError("Please pass a min value, max value, or both")

    def name(self) -> str:
        """Generate a unique name based on all attributes that define this cut."""
        params = [f"column={self.column}"]
        if self.flag is not None:
            params.append(f"flag={hex(self.flag)}")
        if self.min_value is not None:
            params.append(f"min={self.min_value}")
        if self.max_value is not None:
            params.append(f"max={self.max_value}")
        return "Custom Cut (" + ", ".join(params) + ")"

    def __str__(self):
        return self.name()


class UncertaintyCut(Cut):
    def __init__(self, column: str, flag: int = 0x2, max_value: float = 160.0):
        super().__init__(column=column, flag=flag, max_value=max_value)

    @staticmethod
    def name() -> str:
        return "Uncertainty Cut"


class UncertaintyEstimation(Cut):
    def __init__(self, temp_x2_max_value: float = 20, uncert_cut_flag: int = 0x2):
        super().__init__()
        self.temp_x2_max_value = temp_x2_max_value
        self.uncert_cut_flag = uncert_cut_flag

    @staticmethod
    def name() -> str:
        return "Uncertainty Estimation"


class ChiSquareCut(Cut):
    def __init__(
        self,
        column: str,
        flag: int = 0x1,
        max_value: float = 10,
        snr_bound: float = 3,
        min_cut: int = 3,
        max_cut: int = 50,
        cut_step: int = 1,
        use_pre_mjd0_lc: bool = False,
    ):
        super().__init__(column=column, flag=flag, max_value=max_value)

        self.snr_bound = snr_bound

        if min_cut >= max_cut:
            raise ValueError(
                f"Min cut {min_cut} cannot be greater than or equal to max cut {max_cut}"
            )
        if cut_step > (max_cut - min_cut):
            raise ValueError(
                f"Cut step {cut_step} cannot be greater than the difference between max and min cut {max_cut - min_cut}"
            )
        self.min_cut = min_cut
        self.max_cut = max_cut
        self.cut_step = cut_step

        self.use_pre_mjd0_lc = use_pre_mjd0_lc

    @staticmethod
    def name() -> str:
        return "Chi-Square Cut"


class ControlLightCurveCut(Cut):
    def __init__(
        self,
        flag: int = 0x400000,
        questionable_flag: int = 0x80000,
        x2_max: float = 2.5,
        x2_flag: int = 0x100,
        snr_max: float = 3.0,
        snr_flag: int = 0x200,
        Nclip_max: int = 2,
        Nclip_flag: int = 0x400,
        Ngood_min: int = 4,
        Ngood_flag: int = 0x800,
    ):
        super().__init__(flag=flag)

        if flag == questionable_flag:
            raise ValueError(
                f"Bad measurements flag {flag} cannot be equal to questionable measurements flag {questionable_flag}"
            )
        self.questionable_flag = questionable_flag

        self.x2_max = x2_max
        self.x2_flag = x2_flag

        self.snr_max = snr_max
        self.snr_flag = snr_flag

        self.Nclip_max = Nclip_max
        self.Nclip_flag = Nclip_flag

        self.Ngood_min = Ngood_min
        self.Ngood_flag = Ngood_flag

    @staticmethod
    def name() -> str:
        return "Control Light Curve Cut"


class BadDayCut(Cut):
    def __init__(
        self,
        flag: int = 0x800000,
        mjd_bin_size: float = 1.0,
        x2_max: float = 4.0,
        Nclip_max: int = 1,
        Ngood_min: int = 2,
        ixclip_flag: int = 0x1000,
        smallnum_flag: int = 0x2000,
    ):
        super().__init__(flag=flag)
        self.mjd_bin_size = mjd_bin_size
        self.x2_max = x2_max
        self.Nclip_max = Nclip_max
        self.Ngood_min = Ngood_min
        self.ixclip_flag = ixclip_flag
        self.smallnum_flag = smallnum_flag

    @staticmethod
    def name() -> str:
        return "Bad Day Cut"


class CutList:
    def __init__(self):
        self.logger = CustomLogger(self.__class__.__name__)
        self.list: Dict[str, Cut] = {}

    def add(self, cut: Cut):
        if cut.name() in self.list:
            self.logger.warning(
                f"Cut by the name {cut.name()} already exists; overwriting", dots=True
            )
        self.list[cut.name()] = cut

    def add_many(self, cuts: List[Cut]):
        for cut in cuts:
            self.add(cut)

    def get(self, name: str) -> Cut | None:
        if not name in self.list:
            return None
        return self.list[name]

    def remove(self, name: str):
        if self.has(name):
            del self.list[name]

    def remove_many(self, names: List[str]):
        for name in names:
            self.remove(name)

    def remove_by_flag(self, flag: int):
        """
        Removes any Cut object from self.list that has a matching Cut.flag value.
        """
        to_remove = [
            name
            for name, cut in self.list.items()
            if not isinstance(cut, UncertaintyEstimation)
            and cut.flag is not None
            and cut.flag == flag
        ]
        for name in to_remove:
            del self.list[name]

    def has(self, name: str):
        return name in self.list

    def can_apply_directly(self, name: str):
        return self.list[name].can_apply_directly()

    def get_flag_duplicates(self) -> List[int]:
        if len(self.list) < 1:
            return

        unique_flags = set()
        duplicate_flags = []

        for cut in self.list.values():
            if isinstance(cut, UncertaintyEstimation):
                continue

            flags = cut.get_flags()

            for flag in flags:
                if flag in unique_flags:
                    duplicate_flags.append(flag)
                else:
                    unique_flags.add(flag)

        return duplicate_flags

    def get_custom_cuts(self) -> Dict[str, CustomCut]:
        custom_cuts = {}
        for name, cut in self.list.items():
            if isinstance(cut, CustomCut):
                custom_cuts[name] = cut
        return custom_cuts

    def get_all_flags(self):
        mask = 0
        for cut in self.list.values():
            if not isinstance(cut, UncertaintyEstimation):
                flags = cut.get_flags()
                if len(flags) > 0:
                    combined_flags = combine_flags(flags)
                    mask = mask | combined_flags
        return mask

    def get_all_default_flags(self):
        mask = 0
        for cut in self.list.values():
            if not isinstance(cut, UncertaintyEstimation) and cut.flag is not None:
                mask = mask | cut.flag
        return mask

    def get_previous_flags(self, current_cut_name: str):
        skip_names: List = (
            [
                UncertaintyEstimation.name(),
                BadDayCut.name(),
            ]
            + list(self.get_custom_cuts().keys())
            + [
                ControlLightCurveCut.name(),
                ChiSquareCut.name(),
                UncertaintyCut.name(),
            ]
        )

        try:
            current_cut_index = skip_names.index(current_cut_name)
        except:
            raise ValueError(f"No cut by name {current_cut_name} found")
        skip_names = skip_names[: current_cut_index + 1]
        mask = 0
        for name in self.list:
            flag = self.list[name].flag
            if not name in skip_names and flag is not None:
                mask = mask | flag
        return mask

    def iterator(self):
        names: List = (
            [
                UncertaintyCut.name(),
                UncertaintyEstimation.name(),
                ChiSquareCut.name(),
                ControlLightCurveCut.name(),
            ]
            + list(self.get_custom_cuts().keys())
            + [
                BadDayCut.name(),
            ]
        )

        for name in names:
            if name in self.list:
                yield self.list[name]

    def __str__(self):
        output = []
        for name in self.list:
            output.append("• " + self.list[name].__str__())
        return "\n".join(output)
