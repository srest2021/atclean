#!/usr/bin/env python

from abc import ABC, abstractmethod
from configparser import ConfigParser
import configparser
from functools import reduce
from typing import Callable, Dict, Any, List, Optional, Self, Set, Tuple, Type
import re, json, requests, time, sys, io, os
from astropy import units as u
from astropy.coordinates import Angle
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


def AandB(A, B):
    return np.intersect1d(A, B, assume_unique=False)


def AnotB(A, B):
    return np.setdiff1d(A, B)


def AorB(A, B):
    return np.union1d(A, B)


def not_AandB(A, B):
    return np.setxor1d(A, B)


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


def apparent_to_absolute_mag(values, distance_modulus=29.04, precision=2):
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


# load a JSON config file
def load_json_config(filename: str):
    try:
        print(f"Loading JSON config file at {filename}...")
        with open(filename) as cfg:
            return json.load(cfg)
    except Exception as e:
        raise RuntimeError(f"Could not load JSON config file at {filename}: {str(e)}")


def new_row(t: pd.DataFrame, d: Dict = None):
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
    """Parse value from config file by converting None-like string to None."""
    if value:
        stripped = value.strip().lower()
        if stripped == "none":
            return None
        if stripped == "true":
            return True
        if stripped == "false":
            return False
    return value


def parse_comma_separated_string(string: str | None):
    if string is None:
        return None
    return [item.strip() for item in string.split(",")]


def make_dir_if_not_exists(directory):
    """
    Creates a directory if it does not exist. Handles permission errors and other exceptions.

    :param directory: Path to the directory to create.
    """
    if not os.path.isdir(directory):
        try:
            os.makedirs(directory)
        except PermissionError:
            print(f"Permission denied: Cannot create directory at {directory}")
        except FileExistsError:
            # This can occur if the directory is created between the `isdir` check and `makedirs` call.
            print(f"Directory already exists: {directory}")
        except Exception as e:
            print(f"An error occurred while creating directory {directory}: {str(e)}")


# load a .ini config file
def load_config(filename):
    cfg = configparser.ConfigParser()
    try:
        print(f"\nLoading config file at {filename}...")
        cfg.read(filename)
    except Exception as e:
        raise RuntimeError(f"Could not load config file at {filename}: {str(e)}")
    return cfg


def extract_from_subdir(
    directory: str,
    pattern: re.Pattern,
    group_name: str | int,
    convert_function: Callable = lambda x: x,
):
    if not os.path.isdir(directory):
        print(f"WARNING: Cannot search because the path does not exist: {directory}")
        return []

    extracted_values = set()

    for file in os.listdir(directory):
        match = pattern.match(file)
        if match:
            value = match.group(group_name)
            extracted_values.add(convert_function(value))

    if not extracted_values:
        raise RuntimeError(f"Could not find {group_name} from the files in {directory}")

    return extracted_values


def has_match(directory: str, pattern: re.Pattern):
    if not os.path.isdir(directory):
        return False

    for file in os.listdir(directory):
        if pattern.match(file):
            return True
    return False


def find_all_filts(directory: str, tnsname: str):
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


def find_all_control_indices(directory: str, tnsname: str, filt=None):
    if filt is None:
        filt_pattern = r".*"
    else:
        filt_pattern = re.escape(filt)

    pattern = re.compile(
        rf"^{re.escape(tnsname)}_i(?P<index>\d{{3}})"  # captures control index
        rf"\.{filt_pattern}"  # match specific filt if provided
        r"(?:\.\d+\.\d+days)?"  # optional mjdbinsize
        r"(?:\.clean)?"  # optional 'clean'
        r"\.lc\.txt$"  # ends with '.lc.txt'
    )
    subdir = os.path.join(directory, tnsname, "controls")
    return extract_from_subdir(subdir, pattern, "index", convert_function=int)


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


def validate_mjd_ranges(ranges: List[List[int]], var_name: str = "MJD_RANGES") -> None:
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

    print(f"Expanding range edges by {expand_edges}...")

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
    expand_edges: Optional[float] = 0.0,
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
        print(f"Excluding additional MJD ranges {exclude_mjd_ranges}...")
        if len(exclude_mjd_ranges) > 0:
            combined = inverse + exclude_mjd_ranges
            inverse = _merge_ranges(combined)

    return _expand_ranges(inverse, min_mjd, max_mjd, expand_edges=expand_edges)


class PlotLimits:
    def __init__(self, xlower=None, xupper=None, ylower=None, yupper=None):
        self.xlower = xlower
        self.xupper = xupper
        self.ylower = ylower
        self.yupper = yupper

    def set_lims(
        self,
        xlims: Optional[Tuple[float, float]] = None,
        ylims: Optional[Tuple[float, float]] = None,
    ):
        if xlims is not None:
            self.set_xlims(xlims)
        if ylims is not None:
            self.set_ylims(ylims)

    def set_xlims(self, xlims: Tuple[float, float]):
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

    def set_ylims(self, ylims: Tuple[float, float]):
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

    def get_xlims(self):
        if self.xlower is None and self.xupper is None:
            return None
        return self.xlower, self.xupper

    def get_ylims(self):
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
            config_preset_settings: Dict = config[f"column_name_preset.{self.preset}"]
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

        self.optional_columns: Dict[str, Optional[str]] = self._validate_columns_dict(
            {
                "chisquare": config_preset_settings.get("chisquare_column_name"),
                "filt": config_preset_settings.get("filter_column_name"),
                "mag": config_preset_settings.get("mag_column_name"),
                "dmag": config_preset_settings.get("dmag_column_name"),
                "ra": config_preset_settings.get("ra_column_name"),
                "dec": config_preset_settings.get("dec_column_name"),
            }
        )

        # extra columns to copy
        extra_columns = parse_config_str(config_preset_settings["extra_columns"])
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

    def __getattr__(self, name: str) -> Optional[str]:
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

        lines = ["--- Required Columns ---"]
        for k, v in self.required_columns.items():
            if k in skip_colnames:
                continue
            lines.append(f"{k}: {v}")

        lines.append("--- Optional Columns ---")
        for k, v in self.optional_columns.items():
            lines.append(f"{k}: {v}")

        lines.append("--- Extra Columns to Copy ---")
        lines.append(", ".join(self.extra_columns) if self.extra_columns else "(None)")

        return "\n".join(lines)


class Credentials:
    def __init__(
        self, atlas_username, atlas_password, tns_api_key, tns_id, tns_bot_name
    ):
        self.atlas_username = parse_config_str(atlas_username)
        self.atlas_password = parse_config_str(atlas_password)
        self.tns_api_key = parse_config_str(tns_api_key)
        self.tns_id = parse_config_str(tns_id)
        self.tns_bot_name = parse_config_str(tns_bot_name)

    def validate_tns_credentials(self):
        tns_params = [self.tns_api_key, self.tns_id, self.tns_bot_name]
        not_none_count = sum(param is not None for param in tns_params)
        if 0 < not_none_count < 3:
            raise RuntimeError(
                "Either all or none of 'tns_api_key', 'tns_id', and 'tns_bot_name' must be provided."
            )


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
        s = re.compile("\:")
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

    def set_RA(self, ra):
        self.ra = RA(ra)

    def set_Dec(self, dec):
        self.dec = Dec(dec)

    def RA_str(self):
        if self._is_angle_missing(self.ra):
            return np.nan
        return f"{self.ra.angle.degree:0.14f}"

    def Dec_str(self):
        if self._is_angle_missing(self.dec):
            return np.nan
        return f"{self.dec.angle.degree:0.14f}"

    def _is_angle_missing(self, angle: BaseAngle) -> bool:
        return angle.angle is None

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

    def __str__(self):
        output = []
        if self.is_ra_present():
            output.append(f"RA {self.ra.angle.degree:0.14f}")
        if self.is_dec_present():
            output.append(f"Dec {self.dec.angle.degree:0.14f}")

        if len(output) < 1:
            return f"WARNING: Coordinates are empty and cannot be printed."
        return ", ".join(output)


# input/output table containing TNS names, RA, Dec, and MJD0
# (TODO: if MJD0=None, consider entire light curve as pre-SN light curve)
class SnInfoTable:
    def __init__(self, directory, filename=None):
        if filename is None:
            self.filename = f"{directory}/sninfo.txt"
        else:
            self.filename = f"{directory}/{filename}"

        try:
            print(f"Loading SN info table at {self.filename}...")
            self.t = pd.read_table(self.filename, sep="\s+")
            if not "tnsname" in self.t.columns:
                raise RuntimeError('SN info table must have a "tnsname" column.')
            self.t["ra"] = self.t["ra"].astype(str)
            self.t["dec"] = self.t["dec"].astype(str)
            print("Success")
        except Exception:
            print(f"No existing SN info table at that path; creating blank table...")
            self.t = pd.DataFrame(
                columns=["tnsname", "ra", "dec", "mjd0"]
            )  # , 'closebright_ra', 'closebright_dec'])

    def get_row(self, tnsname):
        if self.t.empty:
            # raise RuntimeError(f'Error: Cannot get info for SN {tnsname}--table is empty.')
            return -1, None

        matching_ix = np.where(self.t["tnsname"].eq(tnsname))[0]
        if len(matching_ix) >= 2:
            print(
                f"WARNING: SN info table has {len(matching_ix)} matching rows for TNS name {tnsname}. Dropping duplicate rows..."
            )
            self.t.drop(matching_ix[1:], inplace=True)
            return matching_ix[0], self.t.loc[matching_ix[0], :]
        elif len(matching_ix) == 1:
            return matching_ix[0], self.t.loc[matching_ix[0], :]
        else:
            # raise RuntimeError(f'Error: Cannot get info for SN {tnsname}--row doesn\'t exist.')
            return -1, None

    def is_nan(self, value) -> bool:
        if value is None:
            return True
        if isinstance(value, (float, np.floating)) and np.isnan(value):
            return True
        if isinstance(value, str) and value.strip().lower() in ["nan", "", "none"]:
            return True
        return False

    def get_info(self, tnsname):
        _, row = self.get_row(tnsname)
        if row is None:
            return None, None, None

        # coords = Coordinates(row['ra'], row['dec'])
        ra = None if self.is_nan(row["ra"]) else row["ra"]
        dec = None if self.is_nan(row["dec"]) else row["dec"]

        if self.is_nan(row["mjd0"]):
            mjd0 = None
        else:
            if not isinstance(row["mjd0"], (int, float, np.integer, np.floating)):
                raise RuntimeError(f"Invalid MJD0: {row['mjd0']}")
            mjd0 = float(row["mjd0"])

        return ra, dec, mjd0

    def update_row_at_index(
        self, index, coords: Coordinates = None, mjd0: float = None, overwrite=False
    ):
        try:
            if overwrite or np.isnan(self.t.loc[index, "mjd0"]):
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

    def add_new_row(self, tnsname, coords: Coordinates = None, mjd0: float = None):
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
        self, tnsname, coords: Coordinates = None, mjd0: float = None, overwrite=False
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
            self.add_new_row(tnsname, coords, mjd0)

    def save(self):
        print(f"\nSaving SN info table at {self.filename}...")
        self.t["ra"] = self.t["ra"].astype(str)
        self.t["dec"] = self.t["dec"].astype(str)
        self.t.to_string(self.filename, index=False)
        print("Success")

    def __str__(self):
        return self.t.to_string()


def format_float(value: int | float):
    """
    Format with up to 5 decimals, strip trailing zeros, then ensure at least 1 decimal
    """
    formatted = f"{value:.5f}".rstrip("0").rstrip(".")
    if "." not in formatted:
        formatted += ".0"
    return formatted


def get_filename(
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
        filename += f".{format_float(mjdbinsize)}days"

    if cleaned:
        filename += f".clean"

    filename += ".lc.txt"
    return filename


def query_tns(tnsname, api_key, tns_id, bot_name):
    if tns_id is None or bot_name is None:
        print(
            "WARNING: Cannot query TNS without TNS ID and bot name. Please specify these parameters in config.ini."
        )
        return None

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
        print(json_data["data"])
        raise RuntimeError("ERROR in query_tns(): " + str(e))


def get_tns_coords_from_json(json_data):
    try:
        coords = Coordinates(json_data["data"]["ra"], json_data["data"]["dec"])
        return coords
    except Exception as e:
        raise RuntimeError(f"Failed to get coordinates from TNS JSON data: {str(e)}")


def get_tns_mjd0_from_json(json_data):
    try:
        disc_date = json_data["data"]["discoverydate"]
        date = list(disc_date.partition(" "))[0]
        time = list(disc_date.partition(" "))[2]
        date_object = Time(date + "T" + time, format="isot", scale="utc")
        mjd0 = date_object.mjd - DISC_DATE_BUFFER
        return mjd0
    except Exception as e:
        raise RuntimeError(f"Failed to get discovery date from TNS JSON data: {str(e)}")


def get_mjd0_from_tns(
    tnsname: str, sninfo: SnInfoTable, credentials: Credentials
) -> Tuple[float, Coordinates | None]:
    _, sninfo_row = sninfo.get_row(tnsname)
    if not sninfo_row is None and not np.isnan(sninfo_row["mjd0"]):
        # get MJD0 from SN info table
        print(f'\nSetting MJD0 to {sninfo_row["mjd0"]} MJD from SN info table...')
        mjd0 = float(sninfo_row["mjd0"])
        if not isinstance(mjd0, (int, float)):
            raise RuntimeError(f"Invalid MJD0: {mjd0}")
        else:
            print("Success")
            return mjd0, None
    else:
        # get MJD0 from TNS
        print(f"\nQuerying TNS for SN {tnsname} discovery date...")
        credentials.validate_tns_credentials()
        json_data = query_tns(
            tnsname,
            credentials.tns_api_key,
            credentials.tns_id,
            credentials.tns_bot_name,
        )
        mjd0 = get_tns_mjd0_from_json(json_data)
        coords = get_tns_coords_from_json(json_data)
        return mjd0, coords


def query_atlas(headers, ra, dec, min_mjd, max_mjd):
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
                print(f"Task url: {task_url}")
            elif resp.status_code == 429:
                message = resp.json()["detail"]
                print(f"{resp.status_code} {message}")
                t_sec = re.findall(r"available in (\d+) seconds", message)
                t_min = re.findall(r"available in (\d+) minutes", message)
                if t_sec:
                    waittime = int(t_sec[0])
                elif t_min:
                    waittime = int(t_min[0]) * 60
                else:
                    waittime = 10
                print(f"Waiting {waittime} seconds")
                time.sleep(waittime)
            else:
                print(f"ERROR {resp.status_code}")
                print(resp.text)
                sys.exit()

    result_url = None
    taskstarted_printed = False

    print("Waiting for job to start...")
    while not result_url:
        with requests.Session() as s:
            resp = s.get(task_url, headers=headers)
            if resp.status_code == 200:
                if not (resp.json()["finishtimestamp"] is None):
                    result_url = resp.json()["result_url"]
                    print(f"Task is complete with results available at {result_url}")
                    break
                elif resp.json()["starttimestamp"]:
                    if not taskstarted_printed:
                        print(
                            f"Task is running (started at {resp.json()['starttimestamp']})"
                        )
                        taskstarted_printed = True
                    time.sleep(2)
                else:
                    # print(f"Waiting for job to start (queued at {resp.json()['timestamp']})")
                    time.sleep(4)
            else:
                print(f"ERROR {resp.status_code}")
                print(resp.text)
                sys.exit()

    with requests.Session() as s:
        if result_url is None:
            print("WARNING: Empty light curve (no data within this MJD range).")
            dfresult = pd.DataFrame(
                columns=[
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
            )
        else:
            result = s.get(result_url, headers=headers).text
            dfresult = pd.read_csv(
                io.StringIO(result.replace("###", "")), delim_whitespace=True
            )

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
    print("\nSearching config file for custom cuts...")

    required_keys = {"column", "flag", "max_value", "min_value"}
    custom_cuts = []

    for key in config:
        if key.endswith("_cut") and not key in CONFIG_CUT_NAMES:
            if not required_keys.issubset(config[key].keys()):
                print(
                    f"WARNING: Custom cut {key} missing required fields (required fields: {required_keys}); skipping..."
                )
            else:
                custom_cuts.append(config[key])

    print(f"Found {len(custom_cuts)}")
    return custom_cuts


def get_config_flags(config: ConfigParser) -> List[int]:
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
        column: str = None,
        flag: int = None,
        min_value: float = None,
        max_value: float = None,
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
            and self.column
            and (self.min_value is not None or self.max_value is not None)
        )

    @staticmethod
    @abstractmethod
    def name() -> str:
        pass

    def get_flags(self, keys=False) -> Dict | List:
        """Extract all attributes that end in '_flag' or are 'flag'."""
        flags = {
            key: value
            for key, value in vars(self).items()
            if isinstance(value, int) and (key.endswith("_flag") or key == "flag")
        }
        return flags if keys else list(flags.values())

    def __str__(self) -> str:
        details = []
        if self.column:
            details.append(f"column={self.column}")
        if self.min_value is not None:
            details.append(f"min_value={self.min_value}")
        if self.max_value is not None:
            details.append(f"max_value={self.max_value}")

        flags = self.get_flags(keys=True)
        if flags:
            for flag_name, flag_value in flags.items():
                details.append(f"{flag_name}={hex(flag_value)}")

        return (f"{self.name()}: " + ", ".join(details)) if details else self.name()


class CustomCut(Cut):
    def __init__(
        self, column: str, flag: int, min_value: float = None, max_value: float = None
    ):
        super().__init__(
            column=column, flag=flag, min_value=min_value, max_value=max_value
        )
        if not self.can_apply_directly():
            raise ValueError("Please pass a min value, max value, or both")

    def name(self) -> str:
        """Generate a unique name based on all attributes that define this cut."""
        params = [f"column={self.column}", f"flag={hex(self.flag)}"]
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


class UncertaintyEstimation:
    def __init__(self, temp_x2_max_value: int = 20, uncert_cut_flag: int = 0x2):
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
        self.list: Dict[str, Type[Cut]] = {}

    def add(self, cut: Cut):
        if cut.name() in self.list:
            print(
                f"WARNING: cut by the name {cut.name()} already exists; overwriting..."
            )
        self.list[cut.name()] = cut

    def get(self, name: str):
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

    def get_flag_duplicates(self):
        if len(self.list) < 1:
            return

        unique_flags = set()
        duplicate_flags = []

        for name, cut in self.list.items():
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
            if not isinstance(cut, UncertaintyEstimation):
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
            if not name in skip_names:
                mask = mask | self.list[name].flag
        return mask

    def __str__(self):
        output = ""
        for name in self.list:
            output += self.list[name].__str__()
        return output
