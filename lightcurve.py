#!/usr/bin/env python

from abc import ABC, abstractmethod
from configparser import ConfigParser
from typing import Dict, Any, List, Optional, Self, Set, Tuple, Type
import re, json, requests, time, sys, io
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

DEFAULT_CUT_NAMES = ["uncert_cut", "x2_cut", "controls_cut", "badday_cut", "averaging"]

# ATLAS template change dates
TEMPLATE_CHANGE_1_MJD = 58417
TEMPLATE_CHANGE_2_MJD = 58882

"""
UTILITY
"""


def AandB(A, B):
    return np.intersect1d(A, B, assume_unique=False)


def AnotB(A, B):
    return np.setdiff1d(A, B)


def AorB(A, B):
    return np.union1d(A, B)


def not_AandB(A, B):
    return np.setxor1d(A, B)


def get_allowed_presets(config: ConfigParser) -> list[str]:
    """
    Extract all preset names from the config that match 'column_name_preset.<PRESET NAME>'.
    """
    pattern = re.compile(r"^column_name_preset\.(.+)$")
    return [match.group(1) for key in config.keys() if (match := pattern.match(key))]


def parse_config_value(value: str | None):
    """Parse value from config file by converting None-like string to None."""
    if value and value.strip().lower() == "none":
        return None
    return value


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
            name = parse_config_value(name)
            if no_nones and name is None:
                raise RuntimeError(
                    f"ERROR: Column name '{name}' in config preset {self.preset} cannot be None"
                )
            columns_dict[key] = name
        return columns_dict

    def _read_config(self, config: ConfigParser):
        try:
            config_preset_settings: Dict = config[f"column_name_preset.{self.preset}"]
        except:
            raise RuntimeError(
                f"ERROR: Preset '{self.preset}' (field '{f'column_name_preset.{self.preset}'}') not found in config file."
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
        extra_columns = parse_config_value(config_preset_settings["extra_columns"])
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

        if not overwrite and (
            key in self.required_columns or key in self.optional_columns
        ):
            raise RuntimeError(
                f"ERROR: Column key '{key}' is already defined as {'a required' if key in self.required_columns else 'an optional'} column."
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
                    f"ERROR: Cannot update non-existing required column name {key} with '{name}'"
                )
            self.required_columns[key] = name
        else:
            if key not in self.optional_columns:
                raise RuntimeError(
                    f"ERROR: Cannot update non-existing optional column name {key} with '{name}'"
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
        self.atlas_username = parse_config_value(atlas_username)
        self.atlas_password = parse_config_value(atlas_password)
        self.tns_api_key = parse_config_value(tns_api_key)
        self.tns_id = parse_config_value(tns_id)
        self.tns_bot_name = parse_config_value(tns_bot_name)

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

    def get_RA_str(self):
        if self._is_angle_missing(self.ra):
            return np.nan
        return f"{self.ra.angle.degree:0.14f}"

    def get_Dec_str(self):
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
        filename += f".{mjdbinsize:0.2f}days"

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
        raise RuntimeError(
            f"ERROR: Failed to get coordinates from TNS JSON data: {str(e)}"
        )


def get_tns_mjd0_from_json(json_data):
    try:
        disc_date = json_data["data"]["discoverydate"]
        date = list(disc_date.partition(" "))[0]
        time = list(disc_date.partition(" "))[2]
        date_object = Time(date + "T" + time, format="isot", scale="utc")
        mjd0 = date_object.mjd - DISC_DATE_BUFFER
        return mjd0
    except Exception as e:
        raise RuntimeError(
            f"ERROR: Failed to get discovery date from TNS JSON data: {str(e)}"
        )


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
            self.t = pd.read_table(self.filename, delim_whitespace=True)
            if not "tnsname" in self.t.columns:
                raise RuntimeError('ERROR: SN info table must have a "tnsname" column.')
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

    def is_nan(self, string: str):
        return string.lower() == "nan"

    def get_info(self, tnsname):
        _, row = self.get_row(tnsname)
        if row is None:
            return None, None, None

        # coords = Coordinates(row['ra'], row['dec'])
        ra = None if self.is_nan(row["ra"]) else row["ra"]
        dec = None if self.is_nan(row["dec"]) else row["dec"]

        if np.isnan(row["mjd0"]):
            mjd0 = None
        else:
            if not isinstance(row["mjd0"], (int, float)):
                raise RuntimeError(f'ERROR: Invalid MJD0: {row["mjd0"]}')
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
                f"ERROR: Could not update SN info table at index {index}: {str(e)}"
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
        self.t = pd.concat([self.t, pd.DataFrame([row])], ignore_index=True)

    def update_row(
        self, tnsname, coords: Coordinates = None, mjd0: float = None, overwrite=False
    ):
        if self.t.empty:
            self.add_new_row(tnsname, coords, mjd0)
            return

        matching_ix = np.where(self.t["tnsname"].eq(tnsname))[0]
        if len(matching_ix) > 1:
            raise RuntimeError(
                f"ERROR: SN info table has {len(matching_ix)} matching rows for TNS name {tnsname}."
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


def get_mjd0_from_tns(
    tnsname: str, sninfo: SnInfoTable, credentials: Credentials
) -> Tuple[float, Coordinates | None]:
    _, sninfo_row = sninfo.get_row(tnsname)
    if not sninfo_row is None and not np.isnan(sninfo_row["mjd0"]):
        # get MJD0 from SN info table
        print(f'\nSetting MJD0 to {sninfo_row["mjd0"]} MJD from SN info table...')
        mjd0 = float(sninfo_row["mjd0"])
        if not isinstance(mjd0, (int, float)):
            raise RuntimeError(f"ERROR: Invalid MJD0: {mjd0}")
        else:
            print("Success")
            return mjd0, None
    else:
        # get MJD0 from TNS
        print(f"\nQuerying TNS for SN {tnsname} discovery date...")
        json_data = query_tns(
            tnsname,
            credentials.tns_api_key,
            credentials.tns_id,
            credentials.tns_bot_name,
        )
        mjd0 = get_tns_mjd0_from_json(json_data)
        coords = get_tns_coords_from_json(json_data)
        return mjd0, coords


# TODO (COLUMN): not sure if internal column string should be actual or preset
class Cut:
    def __init__(
        self,
        column: str = None,
        min_value: float = None,
        max_value: float = None,
        flag: int = None,
        params: Dict[str, Any] = None,
    ):
        self.column = column
        self.min_value = min_value
        self.max_value = max_value
        self.flag = flag
        self.params = params

    def can_apply_directly(self):
        if (
            not self.flag
            or not self.column
            or (not self.min_value and not self.max_value)
        ):
            return False
        return True

    def __str__(self):
        output = ""
        if self.column:
            output += f"column={self.column} "
        if self.flag:
            output += f"flag={hex(self.flag)} "
        if self.min_value:
            output += f"min_value={self.min_value} "
        if self.max_value:
            output += f"max_value={self.max_value}"
        return output


class CutList:
    def __init__(self):
        self.list: Dict[str, Type[Cut]] = {}

    def add(self, cut: Cut, name: str):
        if name in self.list:
            raise RuntimeError(f"ERROR: cut by the name {name} already exists.")
        self.list[name] = cut

    def get(self, name: str):
        if not name in self.list:
            return None
        return self.list[name]

    def remove(self, names: str | List[str]):
        if isinstance(names, str):
            if self.has(name):
                del self.list[names]
        else:
            for name in names:
                if self.has(name):
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
            if name == "uncert_est":
                continue

            flags = [cut.flag]
            if cut.params:
                for key in cut.params:
                    if key.endswith("_flag"):
                        flags.append(cut.params[key])

            for flag in flags:
                if flag in unique_flags:
                    if not flag is None:
                        duplicate_flags.append(flag)
                elif not flag is None:
                    unique_flags.add(flag)

        return duplicate_flags

    def get_custom_cuts(self) -> Dict[str, Cut]:
        custom_cuts = {}

        for name in self.list:
            if not name in DEFAULT_CUT_NAMES and name != "uncert_est":
                custom_cuts[name] = self.list[name]

        return custom_cuts

    def get_all_flags(self):
        mask = 0
        for name in self.list:
            if name != "uncert_est":
                mask = mask | self.list[name].flag
        return mask

    def get_previous_flags(self, current_cut_name: str):
        skip_names: List = (
            [
                "uncert_est",
                "badday_cut",
            ]
            + list(self.get_custom_cuts().values())
            + [
                "controls_cut",
                "x2_cut",
                "uncert_cut",
            ]
        )

        try:
            current_cut_index = skip_names.index(current_cut_name)
        except Exception as e:
            raise RuntimeError(
                f"ERROR: Cannot get previous flags for a custom cut: {str(e)}"
            )
        skip_names = skip_names[: current_cut_index + 1]
        mask = 0
        for name in self.list:
            if not name in skip_names:
                mask = mask | self.list[name].flag
        return mask

    def __str__(self):
        output = ""
        for name in self.list:
            output += f"\n{name}: " + self.list[name].__str__()
        return output


"""
LIGHT CURVES
"""


class Supernova:
    def __init__(
        self,
        colnames: PresetColumnNames,
        tnsname: str = None,
        ra: str = None,
        dec: str = None,
        mjd0: float = None,
        filt="o",
    ):
        self.colnames_master = colnames

        self.tnsname = tnsname
        self.coords: Coordinates = Coordinates(ra, dec)
        self.mjd0 = mjd0
        self.filt = filt

        self.lcs: Dict[int, LightCurve] = {}

        self.num_controls = 0
        self.all_indices = None
        self.control_indices = None

    def get(self, control_index=0):
        try:
            return self.lcs[control_index].t
        except:
            raise RuntimeError(
                f"ERROR: Cannot get control light curve {control_index}. Num controls set to {self.num_controls} and {len(self.lcs)} lcs in dictionary."
            )

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
            by=[self.colnames_master.mjd], ignore_index=True, inplace=True
        )

        if self.num_controls == 0:
            return

        if verbose:
            print("\nMaking sure SN and control light curve MJDs match up exactly:")

        sn_sorted_mjd = self.lcs[0].t[self.colnames_master.mjd].to_numpy()

        for control_index in self.get_control_indices():
            # sort by MJD
            self.lcs[control_index].t.sort_values(
                by=[self.colnames_master.mjd], ignore_index=True, inplace=True
            )
            control_sorted_mjd = (
                self.lcs[control_index].t[self.colnames_master.mjd].to_numpy()
            )

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
                                self.colnames_master.mjd: mjd,
                                self.colnames_master.mask: 0,
                            }
                        )

                # remove indices of rows in control light curve for which there is no MJD in the SN lc
                if len(only_control_mjd) > 0:
                    ix_to_skip = []
                    for mjd in only_control_mjd:
                        matching_ix = self.lcs[control_index].ix_equal(
                            self.colnames_master.mjd, mjd
                        )
                        if len(matching_ix) != 1:
                            raise RuntimeError(
                                f"ERROR: Couldn't find MJD={mjd} in MJD column, but should be there!"
                            )
                        ix_to_skip.extend(matching_ix)
                    ix = AnotB(self.lcs[control_index].getindices(), ix_to_skip)
                else:
                    ix = self.lcs[control_index].getindices()

                # sort again
                sorted_ix = self.lcs[control_index].ix_sort_by_cols(
                    self.colnames_master.mjd, indices=ix
                )
                self.lcs[control_index].t = self.lcs[control_index].t.loc[sorted_ix]

            self.lcs[control_index].t.reset_index(drop=True, inplace=True)

        print("Success")

    def prep_for_cleaning(self, verbose=False):
        if verbose:
            print(
                'Adding blank "Mask" columns, replacing infs with NaNs, and calculating flux/dflux...'
            )

        for control_index in self.get_all_indices():
            # add blank 'Mask' column
            self.lcs[control_index].t[self.colnames_master.mask] = 0
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
            raise RuntimeError("ERROR: Cannot apply template correction without MJD0")
        return self.lcs[0].apply_template_correction(
            self.mjd0,
            maskval=maskval,
            region1_offset=region1_offset,
            region2_offset=region2_offset,
            region3_offset=region3_offset,
            num_measurements=num_measurements,
        )

    def apply_cut(self, cut: Cut):
        if not cut.can_apply_directly():
            raise RuntimeError(f"ERROR: Cannot directly apply the following cut: {cut}")

        sn_percent_cut = None
        for control_index in self.get_all_indices():
            percent_cut = self.lcs[control_index].apply_cut(
                cut.column, cut.flag, min_value=cut.min_value, max_value=cut.max_value
            )
            if control_index == 0:
                sn_percent_cut = percent_cut

        return sn_percent_cut

    def get_uncert_est_stats(self, cut: Cut):
        def get_sigma_extra(median_dflux, stdev):
            return max(0, np.sqrt(stdev**2 - median_dflux**2))

        stats = pd.DataFrame(
            columns=["control_index", "median_dflux", "stdev", "sigma_extra"]
        )
        stats["control_index"] = self.get_control_indices()
        stats.set_index("control_index", inplace=True)

        for control_index in self.get_control_indices():
            dflux_clean_ix = self.lcs[control_index].ix_unmasked(
                self.colnames_master.mask, maskval=cut.params["uncert_cut_flag"]
            )
            x2_clean_ix = self.lcs[control_index].ix_inrange(
                colnames=[self.colnames_master.chisquare],
                uplim=cut.params["temp_x2_max_value"],
                exclude_uplim=True,
            )
            clean_ix = AandB(dflux_clean_ix, x2_clean_ix)

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
        for control_index in self.get_all_indices():
            self.lcs[control_index].add_noise_to_dflux(sigma_extra)

    def get_all_controls(self):
        controls = [
            deepcopy(self.lcs[control_index].t)
            for control_index in self.lcs
            if control_index > 0
        ]
        all_controls = LightCurve(self.colnames_master)
        all_controls.t = pd.concat(controls, ignore_index=True)
        return all_controls

    def calculate_control_stats(self, previous_flags):
        print("Calculating control light curve statistics...")

        len_mjd = len(self.lcs[0].t[self.colnames_master.mjd])

        # construct arrays for control lc data
        uJy = np.full((self.num_controls, len_mjd), np.nan)
        duJy = np.full((self.num_controls, len_mjd), np.nan)
        Mask = np.full((self.num_controls, len_mjd), 0, dtype=np.int32)

        i = 1
        for control_index in self.get_control_indices():
            if len(self.lcs[control_index].t) != len_mjd or not np.array_equal(
                self.lcs[0].t[self.colnames_master.mjd],
                self.lcs[control_index].t[self.colnames_master.mjd],
            ):
                raise RuntimeError(
                    f"ERROR: SN lc not equal to control lc for control_index {control_index}! Rerun or debug verify_mjds()."
                )
            else:
                uJy[i - 1, :] = self.lcs[control_index].t[self.colnames_master.flux]
                duJy[i - 1, :] = self.lcs[control_index].t[
                    self.lcs[control_index].colnames.dflux_new
                ]
                Mask[i - 1, :] = self.lcs[control_index].t[self.colnames_master.mask]

            i += 1

        c2_param2columnmapping = self.lcs[0].intializecols4statparams(
            prefix="c2_", format4outvals="{:.2f}", skipparams=["converged", "i"]
        )

        for index in range(uJy.shape[-1]):
            pda4MJD = pdastrostatsclass()
            pda4MJD.t[self.colnames_master.flux] = uJy[0:, index]
            pda4MJD.t[self.lcs[0].colnames.dflux_new] = duJy[0:, index]
            pda4MJD.t[self.colnames_master.mask] = np.bitwise_and(
                Mask[0:, index], previous_flags
            )

            pda4MJD.calcaverage_sigmacutloop(
                self.colnames_master.flux,
                noisecol=self.lcs[0].colnames.dflux_new,
                maskcol=self.colnames_master.mask,
                maskval=previous_flags,
                verbose=1,
                Nsigma=3.0,
                median_firstiteration=True,
            )
            self.lcs[0].statresults2table(
                pda4MJD.statparams, c2_param2columnmapping, destindex=index
            )

    def apply_controls_cut(self, cut: Cut, previous_flags: int):
        self.calculate_control_stats(previous_flags)
        self.lcs[0].t["c2_abs_stn"] = (
            self.lcs[0].t["c2_mean"] / self.lcs[0].t["c2_mean_err"]
        )

        # flag SN measurements
        self.lcs[0].flag_by_control_stats(cut)

        # copy over SN's control cut flags to control light curve 'Mask' columns
        flags_arr = np.full(
            self.lcs[0].t[self.colnames_master.mask].shape,
            (
                cut.flag
                | cut.params["questionable_flag"]
                | cut.params["x2_flag"]
                | cut.params["stn_flag"]
                | cut.params["Nclip_flag"]
                | cut.params["Ngood_flag"]
            ),
        )
        flags_to_copy = np.bitwise_and(
            self.lcs[0].t[self.colnames_master.mask], flags_arr
        )
        for control_index in self.get_control_indices():
            self.lcs[control_index].copy_flags(flags_to_copy)

        # self.drop_extra_columns()

        len_ix = len(self.lcs[0].getindices())
        x2_percent_cut = (
            100
            * len(
                self.lcs[0].ix_masked(
                    self.colnames_master.mask, maskval=cut.params["x2_flag"]
                )
            )
            / len_ix
        )
        stn_percent_cut = (
            100
            * len(
                self.lcs[0].ix_masked(
                    self.colnames_master.mask, maskval=cut.params["stn_flag"]
                )
            )
            / len_ix
        )
        Nclip_percent_cut = (
            100
            * len(
                self.lcs[0].ix_masked(
                    self.colnames_master.mask, maskval=cut.params["Nclip_flag"]
                )
            )
            / len_ix
        )
        Ngood_percent_cut = (
            100
            * len(
                self.lcs[0].ix_masked(
                    self.colnames_master.mask, maskval=cut.params["Ngood_flag"]
                )
            )
            / len_ix
        )
        questionable_percent_cut = (
            100
            * len(
                self.lcs[0].ix_masked(
                    self.colnames_master.mask, maskval=cut.params["questionable_flag"]
                )
            )
            / len_ix
        )
        percent_cut = (
            100
            * len(self.lcs[0].ix_masked(self.colnames_master.mask, maskval=cut.flag))
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

    def apply_badday_cut(self, cut: Cut, previous_flags, flux2mag_sigmalimit=3.0):
        mjdbinsize = cut.params["mjd_bin_size"]
        avg_sn = AveragedSupernova(
            self.colnames_master,
            tnsname=self.tnsname,
            mjd0=self.mjd0,
            filt=self.filt,
            mjdbinsize=mjdbinsize,
        )
        avg_sn.num_controls = self.num_controls
        for control_index in self.get_all_indices():
            avg_sn.set_avg_lc(
                self.lcs[control_index].average(
                    cut,
                    previous_flags,
                    mjdbinsize=mjdbinsize,
                    flux2mag_sigmalimit=flux2mag_sigmalimit,
                ),
                control_index=control_index,
            )

        all_flags = (
            previous_flags
            | cut.flag
            | cut.params["ixclip_flag"]
            | cut.params["smallnum_flag"]
        )
        percent_cut = (
            100
            * len(
                avg_sn.avg_lcs[0].ix_masked(
                    self.colnames_master.mask, maskval=all_flags
                )
            )
            / len(avg_sn.avg_lcs[0].t)
        )
        return avg_sn, percent_cut

    def drop_extra_columns(self):
        for control_index in self.get_all_indices():
            self.lcs[control_index].drop_extra_columns()

    def count_files_in_dir(self, path):
        directory_path = Path(path)
        files = [f for f in directory_path.iterdir() if f.is_file()]
        return len(files)

    def load(self, input_dir, control_index=0, cleaned=False):
        self.lcs[control_index] = LightCurve(
            self.colnames_master, control_index=control_index, filt=self.filt
        )
        self.lcs[control_index].load_lc(input_dir, self.tnsname, cleaned=cleaned)

    def load_all(self, input_dir, num_controls=0, cleaned=False):
        self.lcs = {}
        self.num_controls = 0

        print(f"\nLoading SN light curve and {num_controls} control light curves...")

        # load SN light curve
        self.load(input_dir, cleaned=cleaned)

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
            f"Successfully loaded SN light curve and {self.num_controls} control light curves (control indices: {self.get_control_indices()})"
        )

    def get_all_indices(self):
        if not self.all_indices:
            self.all_indices = list(self.lcs.keys())
            self.all_indices.sort()
        return self.all_indices

    def get_control_indices(self):
        if not self.control_indices:
            self.control_indices = list(self.lcs.keys())
            if 0 in self.control_indices:
                self.control_indices.remove(0)
            self.control_indices.sort()
        return self.control_indices

    def save_all(self, output_dir, overwrite=False, cleaned=True):
        print(
            f'\nDropping extra columns and saving {"cleaned " if cleaned else ""}SN light curve and {self.num_controls} {"cleaned " if cleaned else ""}control light curves...'
        )
        for control_index in self.get_all_indices():
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
    ):
        Supernova.__init__(self, colnames, tnsname, ra, dec, mjd0, filt)
        self.mjdbinsize = mjdbinsize

        self.avg_lcs: Dict[int, AveragedLightCurve] = {}

    def set_avg_lc(self, lc, control_index=0):
        self.avg_lcs[control_index] = deepcopy(lc)

    def set_avg_lcs(self, lcs):
        self.avg_lcs = deepcopy(lcs)

    def get_avg(self, control_index: int = 0):
        try:
            return self.avg_lcs[control_index].t
        except:
            raise RuntimeError(
                f"Cannot get averaged control light curve {control_index}. Num controls set to {self.num_controls} and {len(self.avg_lcs)} lcs in dictionary."
            )

    def load(self, input_dir, control_index=0):
        self.avg_lcs[control_index] = AveragedLightCurve(
            control_index=control_index, filt=self.filt, mjdbinsize=self.mjdbinsize
        )
        self.avg_lcs[control_index].load_lc(input_dir, self.tnsname)

    def load_all(self, input_dir, num_controls=0):
        self.avg_lcs = {}
        self.num_controls = 0

        print(
            f"\nLoading averaged SN light curve and {num_controls} averaged control light curves..."
        )

        # load averaged SN light curve
        self.load(input_dir)

        if num_controls > 0:
            # keep iterating over control indices until we successfully load num_controls averaged light curves
            control_index = 1
            while self.num_controls < num_controls:
                try:
                    self.load(input_dir, control_index=control_index)
                    self.num_controls += 1
                except:
                    print(
                        f"Could not load control light curve {control_index}; skipping..."
                    )
                    del self.avg_lcs[control_index]
                control_index += 1

        print(
            f"Successfully loaded averaged SN light curve and {self.num_controls} averaged control light curves"
        )

    def save_all(self, output_dir, overwrite=False):
        print(
            f"\nDropping extra columns and saving averaged SN light curve and {self.num_controls} averaged control light curves..."
        )
        for control_index in self.get_all_indices():
            self.avg_lcs[control_index].drop_extra_columns()
            self.avg_lcs[control_index].save_lc(
                output_dir, self.tnsname, overwrite=overwrite
            )
        print("Success")

    def get_all_indices(self):
        if not self.all_indices:
            self.all_indices = list(self.avg_lcs.keys())
            self.all_indices.sort()
        return self.all_indices

    def get_control_indices(self):
        if not self.control_indices:
            self.control_indices = list(self.avg_lcs.keys())
            if 0 in self.control_indices:
                self.control_indices.remove(0)
            self.control_indices.sort()
        return self.control_indices

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

    def get_good_indices(self, flag: int):
        return self.ix_unmasked(self.colnames.mask, maskval=flag)

    def get_bad_indices(self, flag: int):
        return self.ix_masked(self.colnames.mask, maskval=flag)

    def can_plot(self, ix: List[int], columns: List[str] = None):
        if columns is None:
            columns = [self.colnames.mjd, self.colnames.flux, self.colnames.dflux_new]

        # check that we are plotting at least one row
        # and that the columns to plot are not all NaN values
        return len(ix) > 0 and not self.t.loc[ix, columns].isna().all().all()

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
            print("Calculating flux/dflux...")
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

    def flag_by_control_stats(self, cut: Cut):
        # flag SN measurements according to given bounds
        flag_x2_ix = self.ix_inrange(
            colnames=["c2_X2norm"], lowlim=cut.params["x2_max"], exclude_lowlim=True
        )
        flag_stn_ix = self.ix_inrange(
            colnames=["c2_abs_stn"], lowlim=cut.params["stn_max"], exclude_lowlim=True
        )
        flag_nclip_ix = self.ix_inrange(
            colnames=["c2_Nclip"], lowlim=cut.params["Nclip_max"], exclude_lowlim=True
        )
        flag_ngood_ix = self.ix_inrange(
            colnames=["c2_Ngood"], uplim=cut.params["Ngood_min"], exclude_uplim=True
        )
        self.update_mask_column(cut.params["x2_flag"], flag_x2_ix)
        self.update_mask_column(cut.params["stn_flag"], flag_stn_ix)
        self.update_mask_column(cut.params["Nclip_flag"], flag_nclip_ix)
        self.update_mask_column(cut.params["Ngood_flag"], flag_ngood_ix)

        # update mask column with control light curve cut on any measurements flagged according to given bounds
        zero_Nclip_ix = self.ix_equal("c2_Nclip", 0)
        unmasked_ix = self.ix_unmasked(
            self.colnames.mask,
            maskval=cut.params["x2_flag"]
            | cut.params["stn_flag"]
            | cut.params["Nclip_flag"]
            | cut.params["Ngood_flag"],
        )
        self.update_mask_column(
            cut.params["questionable_flag"], AnotB(unmasked_ix, zero_Nclip_ix)
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
        self, cut: Cut, previous_flags, mjdbinsize=1.0, flux2mag_sigmalimit=3.0
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
                    cut.params["ixclip_flag"],
                    fluxstatparams["ix_clip"],
                    remove_old=False,
                )

            # if small number within this bin, flag measurements
            if len(range_good_ix) < 3:
                self.update_mask_column(
                    cut.params["smallnum_flag"], range_ix, remove_old=False
                )
                avg_lc.update_mask_column(
                    cut.params["smallnum_flag"], [avglc_index], remove_old=False
                )
            # else check sigmacut bounds and flag
            else:
                is_bad = False
                if fluxstatparams["Ngood"] < cut.params["Ngood_min"]:
                    is_bad = True
                if fluxstatparams["Nclip"] > cut.params["Nclip_max"]:
                    is_bad = True
                if (
                    not (fluxstatparams["X2norm"] is None)
                    and fluxstatparams["X2norm"] > cut.params["x2_max"]
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

    def apply_cut(self, column_name, flag, min_value=None, max_value=None):
        if not column_name in self.t.columns:
            raise RuntimeError(
                f"ERROR: No column name '{column_name}' exists in light curve; cannot apply custom cut"
            )

        all_ix = self.getindices()
        if not min_value is None or not max_value is None:
            kept_ix = self.ix_inrange(
                colnames=[column_name], lowlim=min_value, uplim=max_value
            )
        else:
            raise RuntimeError(
                f"ERROR: Cannot apply cut without min value ({min_value}) or max value ({max_value})."
            )
        cut_ix = AnotB(all_ix, kept_ix)

        self.update_mask_column(flag, cut_ix)

        percent_cut = 100 * len(cut_ix) / len(all_ix)
        return percent_cut

    def update_mask_column(self, flag, indices, remove_old=True):
        if remove_old:
            # remove any old flags of the same value
            self.t[self.colnames.mask] = np.bitwise_and(
                self.t[self.colnames.mask].astype(int), ~flag
            )

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
        if self.colnames.fluxoffset in self.t.columns:
            print("Subtracting previous offset from flux column...")
            self.t[self.colnames.flux] -= self.t[self.colnames.fluxoffset]
        print("Setting current flux offset to 0...")
        self.t[self.colnames.fluxoffset] = 0

    def _update_flux_offset_column(self, offset, region_ix):
        print(f"Adding offset {offset:0.2f} to flux offset column...")
        if not self.colnames.fluxoffset in self.t.columns:
            self.t[self.colnames.fluxoffset] = 0
        self.t.loc[region_ix, self.colnames.fluxoffset] += offset

    def _get_region_mean(self, region_ix, maskval=None) -> float:
        indices = region_ix
        if not maskval is None:
            indices = self.ix_unmasked(
                self.colnames.mask, maskval=maskval, indices=region_ix
            )
        return self.get_mean(self.colnames.flux, indices=indices, round_result=True)

    def _add_offset(self, offset, region_ix):
        print(f"Adding offset {offset:0.2f} to flux column...")
        self.t.loc[region_ix, self.colnames.flux] += offset
        self._update_flux_offset_column(offset, region_ix)

    def _calculate_offset(self, ix1, ix2, num_measurements=40, maskval=None):
        """Calculate the mean difference between two sets of measurements."""
        mean1 = self._get_region_mean(ix1[-num_measurements:], maskval=maskval)
        mean2 = self._get_region_mean(ix2[:num_measurements], maskval=maskval)
        return mean2 - mean1

    def _get_region_indices(self) -> Dict[str : List[int]]:
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
        mjd0,
        region_ix_dict: Dict[str : List[int]],
        maskval=None,
        num_measurements=40,
    ) -> Dict[str : Optional[float]]:
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
            "global": -1 * self._get_region_mean(global_ix, maskval=maskval),
        }

    def _apply_offsets(
        self,
        region_ix_dict: Dict[str : List[int]],
        offset_dict: Dict[str : Optional[float]],
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
        region1_offset=None,
        region2_offset=None,
        region3_offset=None,
    ):
        self._clear_flux_offset_column()

        region_ix_dict = self._get_region_indices()
        offset_dict: Dict[str : Optional[float]] = {
            "1": region1_offset,
            "2": region2_offset,
            "3": region3_offset,
        }

        return self._apply_offsets(region_ix_dict, offset_dict)

    def _auto_template_correction(self, mjd0, maskval=None, num_measurements=40):
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
        region1_offset: Optional[float] = None,
        region2_offset: Optional[float] = None,
        region3_offset: Optional[float] = None,
        num_measurements: int = 40,
    ):
        self.colnames.add("fluxoffset", f"{self.colnames.flux}_offset")
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
                raise RuntimeError(f"ERROR: Missing required column: {column_name}")

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


class LimCutsTable:
    def __init__(self, lc: LightCurve, stn_bound, indices=None):
        self.t = None

        self.lc = lc
        if indices is None:
            indices = self.lc.getindices()
        self.indices = indices

        self.good_ix, self.bad_ix = self.get_goodbad_indices(stn_bound)

    def get_goodbad_indices(self, stn_bound):
        if not self.lc.colnames.fdf in self.lc.t.columns:
            self.lc.calculate_fdf_column()

        good_ix = self.lc.ix_inrange(
            colnames=[self.lc.colnames.fdf],
            lowlim=-stn_bound,
            uplim=stn_bound,
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
            self.t = pd.concat([self.t, pd.DataFrame([row])], ignore_index=True)


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
                self.mjd = get_tns_mjd0_from_json(json_data)
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
            raise RuntimeError(
                f"ERROR: max MJD {max_mjd} cannot be than min MJD {min_mjd}."
            )

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
    def save(self, input_dir, tnsname, overwrite=False):
        if self.t is None:
            raise RuntimeError(
                "ERROR: Cannot save light curve that hasn't been downloaded yet."
            )

        lc = LightCurve(control_index=self.control_index)
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
        self.peak_appmag = None

    @abstractmethod
    def get_sim_flux(self, mjds, peak_appmag, **kwargs):
        """
        Compute the simulated flux for the given MJDs and peak apparent magnitude.

        :param mjds: List or array of MJDs.
        :param peak_appmag: Desired peak apparent magnitude of the simulation.

        :return: An array of flux values corresponding to the input MJDs.
        """
        pass

    def __str__(self):
        return f'Simulation with model name "{self.model_name}": peak appmag = {self.peak_appmag:0.2f}'


class SimDetecSupernova(AveragedSupernova):
    def __init__(
        self,
        colnames: PresetColumnNames,
        tnsname: str = None,
        mjdbinsize: float = 1.0,
        filt: str = "o",
    ):
        AveragedSupernova.__init__(
            self, colnames, tnsname=tnsname, mjdbinsize=mjdbinsize, filt=filt
        )
        self.avg_lcs: Dict[int, SimDetecLightCurve] = {}

    def apply_rolling_sums(self, sigma_kern: float, flag=0x800000):
        for control_index in self.get_all_indices():
            self.avg_lcs[control_index].apply_rolling_sum(sigma_kern, flag=flag)

    def remove_rolling_sums(self):
        for control_index in self.get_all_indices():
            self.avg_lcs[control_index].remove_rolling_sum()

    def remove_simulations(self):
        for control_index in self.get_all_indices():
            self.avg_lcs[control_index].remove_simulations()

    def load(self, input_dir, control_index=0):
        self.avg_lcs[control_index] = SimDetecLightCurve(
            self.colnames_master,
            control_index=control_index,
            filt=self.filt,
            mjdbinsize=self.mjdbinsize,
        )
        if control_index == 0:
            filename = f"{input_dir}/{self.tnsname}.{self.filt}.{self.mjdbinsize:0.2f}days.lc.txt"
        else:
            filename = f"{input_dir}/controls/{self.tnsname}_i{control_index:03d}.{self.filt}.{self.mjdbinsize:0.2f}days.lc.txt"
        self.avg_lcs[control_index].load_lc_by_filename(filename)


class SimDetecLightCurve(AveragedLightCurve):
    def __init__(
        self,
        colnames: PresetColumnNames,
        control_index=0,
        filt="o",
        mjd0=None,
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
        self.pre_mjd0_ix = self.ix_inrange(self.colnames.mjd, uplim=mjd0)
        self.valid_seasons_ix = None

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

    # apply a rolling sum to the light curve and add SNR, SNRsum, and SNRsumnorm columns
    def apply_rolling_sum(self, sigma_kern, indices=None, flag=0x800000, verbose=False):
        if indices is None:
            indices = self.getindices()
        if len(indices) < 1:
            raise RuntimeError(
                "ERROR: not enough measurements to apply simulated gaussian"
            )
        good_ix = AandB(indices, self.ix_unmasked(self.colnames.mask, flag))

        self.remove_rolling_sum()
        self.cur_sigma_kern = sigma_kern
        self.t.loc[indices, self.colnames.snr] = 0.0
        self.t.loc[good_ix, self.colnames.snr] = (
            self.t.loc[good_ix, self.colnames.flux]
            / self.t.loc[good_ix, self.colnames.dflux]
        )

        new_gaussian_sigma = round(sigma_kern / self.mjdbinsize)
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

    # add simulated flux to the light curve and add SNRsim and SNRsimsum columns
    def add_sim_flux(
        self,
        good_ix,
        sim_flux,
        cur_sigma_kern=None,
        verbose=False,
        remove_old=True,
    ):
        """
        Add simulated flux to the light curve ("uJysim" column) and add "SNRsim" and "SNRsimsum" columns.

        :param lc: Light curve to add the simulated flux to.
        :param good_ix: Unmasked/unflagged indices of the light curve.
        :param cur_sigma_kern: The current kernel size of the rolling sum.
        :param remove_old: Remove any old simulations before adding the simulated flux.
        """
        if cur_sigma_kern is None:
            cur_sigma_kern = self.cur_sigma_kern
        if cur_sigma_kern is None:
            raise RuntimeError(
                "ERROR: No current sigma kern passed as argument or stored during previously applied rolling sum."
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

    # add any simulation to the light curve, specifying parameters using keyword arguments
    def add_simulation(
        self,
        sim: Simulation,
        peak_appmag: float,
        cur_sigma_kern: int = None,
        flag: int = 0x800000,
        verbose: bool = False,
        remove_old: bool = True,
        **kwargs,
    ) -> Self:
        """
        Add any Simulation object to a copy of the light curve, specifying parameters using keyword arguments.

        :param sim: The Simulation to add.
        :param peak_appmag: The desired peak apparent magnitude of the Simulation to add.
        :param cur_sigma_kern: The current sigma of the rolling sum.
        :param flag: The flag value by which to filter out any flagged bins.
        :param remove_old: Remove any old simulations before adding the simulated flux.
        """
        if verbose:
            print(f"Adding simulation: {sim}")

        lc = deepcopy(self)
        good_ix = AandB(lc.getindices(), lc.ix_unmasked(self.colnames.mask, flag))
        sim_flux = sim.get_sim_flux(
            lc.t.loc[good_ix, self.colnames.mjd], peak_appmag, **kwargs
        )

        return lc.add_sim_flux(
            good_ix,
            sim_flux,
            cur_sigma_kern=cur_sigma_kern,
            verbose=verbose,
            remove_old=remove_old,
        )

    # get max FOM (for simulated FOM, column=SNRsimsum; else column=SNRsumnorm)
    # of measurements within the given indices
    def get_max_fom(self, indices=None):
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
