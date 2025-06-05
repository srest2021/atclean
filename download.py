#!/usr/bin/env python

"""
Code adapted from Qinan Wang and Armin Rest by Sofia Rest

Possible inputs:
- .txt table with TNS names and either (RA, Dec, and MJD0) or (TNS bot credentials)
- list of TNS names in command line and either (.txt table with TNS names, RA, Dec, and MJD0) or (TNS credentials)

Outputs:
- downloaded light curve files
- if using TNS credentials, new or updated .txt table with TNS names, RA, Dec, and MJD0
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Type
import os, sys, requests, argparse, configparser, math
import pandas as pd
import numpy as np
from getpass import getpass
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from lightcurve import FullLightCurve, LightCurve
from utils import (
    DISC_DATE_BUFFER,
    Coordinates,
    Credentials,
    SnInfoTable,
    find_all_control_indices,
    find_all_filts,
    get_tns_data,
    is_sn_in_subdir,
    load_config,
    make_dir_if_not_exists,
    new_row,
    parse_comma_separated_string,
)

CTRL_COORDINATES_COLNAMES = [
    "tnsname",
    "control_index",
    "ra",
    "dec",
    "ra_offset",
    "dec_offset",
    "radius_arcsec",
    "n_detec",
    "n_detec_c",
    "n_detec_o",
]


class ControlCoordinatesTable:
    def __init__(self):
        self.t: Optional[pd.DataFrame] = None
        self.num_controls: Optional[int] = None
        self.radius: Optional[Angle] = None
        self.sn_coords: Optional[Coordinates] = None
        self.center_coords: Optional[Coordinates] = None
        self._closebright: bool = False
        self.sn_min_dist: Optional[float] = None

    def init_load(self, directory: str, tnsname: str):
        filename = self.get_filepath(directory, tnsname)
        print(f"Loading control coordinates table at {filename}...")
        self.load(filename)
        print("✅ Success")
        print(self.__str__())

        self._set_num_controls_from_t()

    def init_read_from_file(self, filepath: str):
        print(f"Loading control coordinates table at {filepath}...")
        self._read(filepath)
        print("✅ Success")
        print(self.__str__())

        self._set_num_controls_from_t()

    def init_closebright(
        self, center_coords: Coordinates, num_controls: int, sn_min_dist: float
    ):
        self.num_controls = num_controls
        # set self.radius later

        self._closebright = True
        self.center_coords = center_coords
        self.sn_min_dist = sn_min_dist

        print(
            f"Setting circle pattern of {self.num_controls} control light curves around center location {self.center_coords},"
            f' with minimum {self.sn_min_dist}" distance from SN'
        )

    def init_default(self, num_controls: int, radius: float):
        self.num_controls = num_controls
        self.radius = Angle(radius, u.arcsec)
        # set self.center_coords later

        print(
            f'Setting circle pattern of {self.num_controls} control light curves around SN location with radius of {self.radius}" from center'
        )

    def ready_for_construction(self):
        assert self.num_controls is not None
        assert self.radius is not None
        assert self.sn_coords is not None
        assert self.center_coords is not None

        if self._closebright:
            assert self.sn_min_dist is not None

    def reset_table(self):
        self.t = pd.DataFrame(columns=CTRL_COORDINATES_COLNAMES)

    def _set_num_controls_from_t(self):
        if self.t is None:
            raise RuntimeError("Table (self.t) cannot be None")
        # first row is SN, remaining are controls
        self.num_controls = len(self.t) - 1

    def _read(self, filepath: str):
        try:
            self.t = pd.read_table(filepath, sep="\s+")
            if not "ra" in self.t.columns or not "dec" in self.t.columns:
                raise RuntimeError(
                    'Control coordinates table must have "ra" and "dec" columns.'
                )
        except Exception as e:
            raise RuntimeError(
                f"Could not load control coordinates table at {filepath}: {str(e)}"
            )

        self.num_controls = len(self.t)
        if (
            "control_index" not in self.t.columns
            or self.t["control_index"].isnull().all()
        ):
            self.t["control_index"] = range(1, self.num_controls + 1)

        for colname in CTRL_COORDINATES_COLNAMES:
            if not colname in self.t.columns:
                self.t[colname] = np.full(len(self.t), np.nan)

    def load(self, directory, tnsname, filename: Optional[str] = None):
        if filename is None:
            if tnsname is None:
                raise RuntimeError(
                    "Please provide either a filename or a TNS name to save the control coordinates table."
                )
            filepath = self.get_filepath(directory, tnsname)
        else:
            filepath = f"{directory}/{tnsname}/{filename}"

        if not os.path.exists(filepath):
            raise ValueError(f"File path {filepath} does not exist")

        try:
            self.t = pd.read_table(filepath, sep="\s+")
            if not "ra" in self.t.columns or not "dec" in self.t.columns:
                raise RuntimeError(
                    'Control coordinates table must have "ra" and "dec" columns.'
                )
        except Exception as e:
            raise RuntimeError(
                f"Could not load control coordinates table at {filepath}: {str(e)}"
            )

    def update_filt_lens(
        self, control_index: int, total_len: int, filt_lens: Dict[str, int]
    ):
        if self.t is None:
            raise RuntimeError("Table (self.t) cannot be None")

        indices = np.where(self.t["control_index"] == control_index)[0]
        if len(indices) == 0:
            raise RuntimeError(
                f"Cannot update row in control coordinates table for control index {control_index}: no matching rows."
            )
        elif len(indices) > 1:
            raise RuntimeError(
                f"Cannot update row in control coordinates table for control index {control_index}: duplicate rows."
            )
        index = indices[0]

        # update corresponding row in table with total and filter counts
        # total_len, filt_lens = full_control_lc.get_filt_lens()
        self.t.at[index, "n_detec"] = total_len
        for filt in filt_lens:
            self.t.at[index, f"n_detec_{filt}"] = filt_lens[filt]

    def add_row(
        self,
        tnsname: str,
        control_index: int,
        coords: Coordinates,
        ra_offset: float | Angle = 0.0,
        dec_offset: float | Angle = 0.0,
        radius: float | Angle = 0.0,
        n_detec: int = 0,
        filt_lens: Optional[Dict[str, int]] = None,
    ):
        row = {
            "tnsname": tnsname,
            "control_index": control_index,
            "ra": coords.get_RA_str(),
            "dec": coords.get_Dec_str(),
            "ra_offset": (
                f"{ra_offset.degree:0.14f}"
                if isinstance(ra_offset, Angle)
                else ra_offset
            ),
            "dec_offset": (
                f"{dec_offset.degree:0.14f}"
                if isinstance(dec_offset, Angle)
                else dec_offset
            ),
            "radius_arcsec": radius.arcsecond if isinstance(radius, Angle) else radius,
            "n_detec": n_detec,
        }

        if not filt_lens is None:
            for filt in filt_lens:
                row[f"n_detec_{filt}"] = filt_lens[filt]

        self.t = new_row(self.t, row)

    def construct_row(self, control_index: int):
        if self.num_controls is None:
            raise RuntimeError("Number of control light curves cannot be None")

        angle = Angle(control_index * 360.0 / self.num_controls, u.degree)

        ra_distance = Angle(self.radius.degree * math.cos(angle.radian), u.degree)
        ra_offset = Angle(
            ra_distance.degree * (1.0 / math.cos(self.center_coords.dec.angle.radian)),
            u.degree,
        )
        ra = Angle(self.center_coords.ra.angle.degree + ra_offset.degree, u.degree)

        dec_offset = Angle(self.radius.degree * math.sin(angle.radian), u.degree)
        dec = Angle(self.center_coords.dec.angle.degree + dec_offset.degree, u.degree)

        coords = Coordinates()
        coords.ra.angle = ra
        coords.dec.angle = dec

        if self._closebright:
            # check to see if control light curve location is within minimum distance from SN location
            offset_sep = self.sn_coords.get_distance(coords).arcsecond
            if offset_sep < self.sn_min_dist:
                print(
                    f'Control light curve {control_index:3d} too close to SN location ({offset_sep}" away) with minimum distance to SN as {self.sn_min_dist}"; skipping control light curve...'
                )
                return

        self.add_row(
            str(np.nan),
            control_index,
            coords,
            ra_offset=ra_offset,
            dec_offset=dec_offset,
            radius=self.radius,
        )

    def construct(
        self,
        tnsname: str,
        sn_coords: Coordinates,
    ):
        self.sn_coords = sn_coords
        self.reset_table()

        # add row for SN position
        if self._closebright:
            # circle pattern radius is distance between SN and close bright object
            self.radius = self.sn_coords.get_distance(self.center_coords)

            self.add_row(
                tnsname,
                0,
                self.sn_coords,
                ra_offset=np.nan,
                dec_offset=np.nan,
                radius=self.radius,
            )
        else:
            # center coordinates are the SN location
            self.center_coords = self.sn_coords

            self.add_row(
                tnsname,
                0,
                self.sn_coords,
            )

        self.ready_for_construction()

        # add row for each control light curve
        for i in range(1, self.num_controls + 1):
            self.construct_row(i)

        print("Control light curve coordinates generated: \n", self.__str__())

    def iterator(self, include_sn: bool = False):
        """
        Yields tuples of (control_index, coordinates) from the control coordinates table.

        :param include_sn (bool): If True, include the first row (SN). Defaults to False.
        """
        if self.t is None:
            raise RuntimeError("ControlCoordinatesTable is empty (self.t is None)")

        df = self.t if include_sn else self.t.iloc[1:]
        for _, row in df.iterrows():
            yield row["control_index"], Coordinates(ra=row["ra"], dec=row["dec"])

    def get_filepath(self, directory, tnsname):
        return os.path.join(directory, tnsname, f"{tnsname}_control_coords_table.txt")

    def get_formatters(self):
        def format_int_or_nan(x):
            if pd.isna(x):
                return "NaN"
            return f"{int(x)}"

        int_cols = [col for col in self.t.columns if col.startswith("n_detec")]
        return {col: format_int_or_nan for col in int_cols}

    def save(
        self,
        directory: str,
        tnsname: str,
        filename: Optional[str] = None,
        overwrite: bool = False,
    ):
        if filename is None:
            if tnsname is None:
                raise RuntimeError(
                    "Please provide either a filename or a TNS name to save the control coordinates table."
                )
            filepath = self.get_filepath(directory, tnsname)
        else:
            filepath = f"{directory}/{tnsname}/{filename}"

        print(f"💾 Saving control coordinates table at {filepath}...")
        if self.t is None:
            raise RuntimeError(
                "Cannot save ControlCoordinatesTable: table (self.t) is None"
            )

        if overwrite or not os.path.exists(filepath):
            with pd.option_context("display.float_format", "{:,.14f}".format):
                self.t.to_string(
                    filepath,
                    index=False,
                    formatters=self.get_formatters(),
                )

    def __str__(self):
        if self.t is None:
            return "ControlCoordinatesTable is empty (no data loaded)"
        with pd.option_context("display.float_format", "{:,.8f}".format):
            return self.t.to_string(formatters=self.get_formatters())


class ControlCoordinatesTableFactory:
    @staticmethod
    def validate(
        download_controls: bool,
        num_controls: int = 0,
        control_coords_table_filepath: Optional[str] = None,
        radius: Optional[float] = None,
        center_coords: Optional[Coordinates] = None,
        sn_min_dist: Optional[float] = None,
    ):
        """
        Validates the input arguments for constructing a ControlCoordinatesTable.
        Ensures consistency between download mode and related parameters.
        """
        if not download_controls and (
            control_coords_table_filepath is not None
            or num_controls != 0
            or (center_coords is not None and center_coords.is_complete())
        ):
            raise ValueError(
                "If not downloading control light curves, none of the following should be provided: "
                "File path to table of control light curve coordinates (`control_coords_table_filepath`), nonzero number of control light curves to download (`num_controls`), or center coordinates of the circle pattern (`center_coords`)"
            )

        if download_controls:
            if radius is None or num_controls < 1:
                raise ValueError(
                    "If downloading control light curves in circle pattern, both radius (`radius`) and number of controls to download (`num_controls`) must be provided"
                )

            if (
                center_coords is not None
                and center_coords.is_complete()
                and sn_min_dist is None
            ):
                raise ValueError(
                    "If centering circle pattern around a different location, both coordinates of the new center (`center_coords`) and minimum distance from SN location (`sn_min_dist`) must be provided"
                )

    @staticmethod
    def new(
        download_controls: bool,
        num_controls: int = 0,
        control_coords_table_filepath: Optional[str] = None,
        radius: Optional[float] = None,
        center_coords: Optional[Coordinates] = None,
        sn_min_dist: Optional[float] = None,
    ) -> ControlCoordinatesTable:
        ControlCoordinatesTableFactory.validate(
            download_controls,
            num_controls=num_controls,
            control_coords_table_filepath=control_coords_table_filepath,
            radius=radius,
            center_coords=center_coords,
            sn_min_dist=sn_min_dist,
        )

        control_coords_table = ControlCoordinatesTable()

        print()
        if control_coords_table_filepath:
            control_coords_table.init_read_from_file(control_coords_table_filepath)
        elif (
            center_coords is not None
            and center_coords.is_complete()
            and sn_min_dist is not None
        ):
            control_coords_table.init_closebright(
                center_coords, num_controls, sn_min_dist
            )
        else:
            control_coords_table.init_default(num_controls, radius)

        return control_coords_table


class AtlasAuthenticator:
    @staticmethod
    def authenticate(username: str, password: str) -> Dict[str, str]:
        print("\nConnecting to ATLAS API...")
        resp = requests.post(
            url=f"https://fallingstar-data.com/forcedphot/api-token-auth/",
            data={"username": username, "password": password},
        )
        if resp.status_code != 200:
            raise RuntimeError(f"Authentication failed: {resp.status_code}")
        token = resp.json()["token"]
        print(f"Token: {token}")
        headers = {"Authorization": f"Token {token}", "Accept": "application/json"}
        return headers


def parse_arg_coords(arg_coords: Optional[str]) -> Optional[Coordinates]:
    parsed_coords = parse_comma_separated_string(arg_coords)
    if parsed_coords is None:
        return None
    if len(parsed_coords) > 2:
        raise RuntimeError(
            "Too many coordinates in argument! Please provide comma-separated RA and Dec onlyy."
        )
    if len(parsed_coords) < 2:
        raise RuntimeError(
            "Too few coordinates in argument! Please provide comma-separated RA and Dec."
        )

    return Coordinates(parsed_coords[0], parsed_coords[1])


def resolve_sn_coords_and_mjd0(
    tnsname: str,
    sninfo: Optional[SnInfoTable] = None,
    creds: Optional[Credentials] = None,
    arg_mjd0: Optional[float] = None,
    arg_sn_coords: Optional[Coordinates] = None,
    arg_center_coords: Optional[Coordinates] = None,
    use_disc_date_buffer: bool = True,
) -> tuple[float, Coordinates]:
    print("\n--- Resolving SN coordinates, center coordinates, and MJD0 ---")
    sn_coords, center_coords, mjd0 = None, None, None

    # first, try SN info table
    if sninfo is not None:
        sn_coords, center_coords, mjd0 = sninfo.get_info(tnsname)
        print(
            f"From SnInfoTable:\n"
            f"  SN coordinates: {sn_coords}\n"
            f"  Center coordinates: {center_coords}\n"
            f"  MJD0: {mjd0} MJD"
        )

    # next, overwrite defaults with command line args
    if arg_sn_coords is not None:
        sn_coords = arg_sn_coords
        print(f"Overriding SnInfoTable SN coordinates with --sn_coords: {sn_coords}")
    if arg_center_coords is not None:
        center_coords = arg_center_coords
        print(
            f"Overriding SnInfoTable center coordinates with --center_coords: {center_coords}"
        )
    if arg_mjd0 is not None:
        mjd0 = arg_mjd0
        print(f"Overriding SnInfoTable MJD0 with --mjd0: {mjd0} MJD")

    # now try querying TNS for missing info
    if sn_coords is None or sn_coords.is_incomplete() or mjd0 is None or np.isnan(mjd0):
        if creds is None:
            raise ValueError(
                "Cannot find coordinates or MJD0 in command line or SnInfoTable, but TNS credentials not provided"
            )
        creds.validate_tns_credentials()

        tns_mjd0, tns_sn_coords = get_tns_data(
            tnsname,
            creds.tns_api_key,
            creds.tns_id,
            creds.tns_bot_name,
            use_disc_date_buffer=use_disc_date_buffer,
        )

        if sn_coords is None or sn_coords.is_incomplete():
            sn_coords = tns_sn_coords
            print(f"Using SN coordinates from TNS API: {sn_coords}")

        if mjd0 is None or np.isnan(mjd0):
            mjd0 = tns_mjd0
            print(
                f"Using MJD0 from TNS discovery date{f' - buffer of {DISC_DATE_BUFFER} MJD' if use_disc_date_buffer else ''}: {mjd0} MJD"
            )

    # make sure nothing is missing
    if sn_coords is None or sn_coords.is_incomplete():
        raise RuntimeError(
            "Could not resolve SN coordinates from command line, SnInfoTable, or TNS"
        )
    if mjd0 is None or np.isnan(mjd0):
        raise RuntimeError(
            "Could not resolve SN MJD0 from command line, SnInfoTable, or TNS discovery date"
        )

    return sn_coords, center_coords, mjd0


class AtlasLightCurveDownloader:
    def __init__(
        self, atclean_input_dir: str, atlas_username: str, atlas_password: str
    ):
        self.atclean_input_dir = atclean_input_dir

        self.headers = AtlasAuthenticator.authenticate(atlas_username, atlas_password)
        if self.headers is None:
            raise RuntimeError("No token header!")

        self._lcs: Dict[int, FullLightCurve] = None

    @property
    def lcs(self):
        if self._lcs is None:
            self._lcs = {}
        return self._lcs

    def download_lc(
        self,
        control_index: int,
        coords: Coordinates,
        lookbacktime: Optional[float] = None,
        max_mjd: Optional[float] = None,
    ) -> tuple[int, Dict[str, int]]:
        self.lcs[control_index] = FullLightCurve(
            control_index, ra=coords.get_RA_str(), dec=coords.get_Dec_str()
        )
        self.lcs[control_index].download(
            self.headers, lookbacktime=lookbacktime, max_mjd=max_mjd
        )
        return self.lcs[control_index].get_filt_lens()

    def load_existing_lc(
        self, tnsname: str, control_index: int
    ) -> tuple[int, Dict[str, int]]:
        filt_lens = {}

        for filt in find_all_filts(self.atclean_input_dir, tnsname):
            lc = LightCurve(None, control_index=control_index, filt=filt)
            lc.load_lc(self.atclean_input_dir, tnsname)
            filt_lens[filt] = len(lc.t)

        total_len = sum([filt_len for filt_len in filt_lens.values()])

        return total_len, filt_lens

    def update_and_save_control_coords_table(
        self,
        control_coords_table: ControlCoordinatesTable,
        tnsname: str,
        control_index: int,
        total_len: int,
        filt_lens: Dict[str, int],
    ) -> ControlCoordinatesTable:
        """
        Update and save the current ControlCoordinatesTable,
        so that if it crashes later, we can just load it again.
        """
        control_coords_table.update_filt_lens(control_index, total_len, filt_lens)
        control_coords_table.save(self.atclean_input_dir, tnsname, overwrite=True)
        return control_coords_table

    def download_and_save_sn_lc(
        self,
        tnsname: str,
        sn_coords: Coordinates,
        control_coords_table: ControlCoordinatesTable,
        lookbacktime: Optional[float] = None,
        max_mjd: Optional[float] = None,
        overwrite: bool = False,
    ) -> ControlCoordinatesTable:
        if not overwrite and is_sn_in_subdir(self.atclean_input_dir, tnsname):
            print(
                f"Overwrite set to False and SN light curve already exists; skipping download..."
            )
            if control_coords_table.t is None:
                # load previously saved ControlCoordinatesTable
                control_coords_table = ControlCoordinatesTable()
                control_coords_table.load(input_dir, tnsname)

            print(
                "Control light curve coordinates table loaded: \n",
                control_coords_table,
                "\n--- Download: SN light curve ---",
                "\nSkipped",
            )

            return control_coords_table

        control_coords_table.construct(tnsname, sn_coords)
        print()

        print("\n--- Download: SN light curve ---")
        total_len, filt_lens = self.download_lc(
            0, sn_coords, lookbacktime=lookbacktime, max_mjd=max_mjd
        )
        self.save_downloaded_lc(tnsname, 0, overwrite=overwrite)

        return self.update_and_save_control_coords_table(
            control_coords_table, tnsname, 0, total_len, filt_lens
        )

    def download_and_save_control_lcs(
        self,
        tnsname: str,
        control_coords_table: ControlCoordinatesTable,
        lookbacktime: Optional[float] = None,
        max_mjd: Optional[float] = None,
        overwrite: bool = False,
    ) -> ControlCoordinatesTable:
        existing_control_indices = find_all_control_indices(
            self.atclean_input_dir, tnsname
        )

        for control_index, coords in control_coords_table.iterator():
            print(f"\n--- Download: Control light curve {control_index} ---")
            if not overwrite and control_index in existing_control_indices:
                print(
                    f"Overwrite set to {overwrite} and light curve already exists; skipping download..."
                )
                total_len, filt_lens = self.load_existing_lc(tnsname, control_index)
            else:
                total_len, filt_lens = self.download_lc(
                    control_index, coords, lookbacktime=lookbacktime, max_mjd=max_mjd
                )
                self.save_downloaded_lc(tnsname, control_index, overwrite=overwrite)

            control_coords_table = self.update_and_save_control_coords_table(
                control_coords_table, tnsname, control_index, total_len, filt_lens
            )

        return control_coords_table

    def download_and_save(
        self,
        tnsname: str,
        sn_coords: Coordinates,
        control_coords_table: ControlCoordinatesTable,
        lookbacktime: Optional[float] = None,
        max_mjd: Optional[float] = None,
        overwrite: bool = False,
    ) -> ControlCoordinatesTable:
        control_coords_table = self.download_and_save_sn_lc(
            tnsname,
            sn_coords,
            control_coords_table,
            lookbacktime=lookbacktime,
            max_mjd=max_mjd,
            overwrite=overwrite,
        )
        assert control_coords_table is not None and control_coords_table.t is not None

        control_coords_table = self.download_and_save_control_lcs(
            tnsname,
            control_coords_table,
            lookbacktime=lookbacktime,
            max_mjd=max_mjd,
            overwrite=overwrite,
        )

        return control_coords_table

    def save_downloaded_lc(
        self, tnsname: str, control_index: str, overwrite: bool = False
    ):
        if control_index not in self.lcs.keys():
            raise ValueError(
                f"Light curve with control_index {control_index} has not been downloaded yet"
            )
        self.lcs[control_index].save(
            None, self.atclean_input_dir, tnsname, overwrite=overwrite
        )

    def save_downloaded_lcs(self, tnsname: str, overwrite: bool = False):
        for control_index in self.lcs.keys():
            self.save_downloaded_lc(tnsname, control_index, overwrite=overwrite)


# define command line arguments
def define_args(parser=None, usage=None, conflict_handler="resolve"):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

    parser.add_argument(
        "tnsnames", nargs="+", help="TNS names of the objects to download from ATLAS"
    )
    parser.add_argument(
        "--sninfo_file",
        default=None,
        type=str,
        help="file name of .txt file with SN info table",
    )
    parser.add_argument(
        "--config_file",
        default="config.ini",
        type=str,
        help="file name of .ini file with settings for this class",
    )
    parser.add_argument(
        "-l", "--lookbacktime", default=None, type=int, help="lookback time (MJD)"
    )
    # parser.add_argument('--min_mjd', default=None, type=float, help='minimum MJD to download')
    parser.add_argument(
        "--max_mjd", default=None, type=float, help="maximum MJD to download"
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        default=False,
        action="store_true",
        help="overwrite existing file with same file name",
    )

    # for downloading single SN only
    parser.add_argument(
        "--sn_coords",
        type=parse_arg_coords,
        default=None,
        help="comma-separated RA and Dec of SN light curve to download",
    )
    parser.add_argument(
        "--mjd0", type=float, default=None, help="transient start date in MJD"
    )

    # for control light curves
    parser.add_argument(
        "-c",
        "--download_controls",
        default=False,
        action="store_true",
        help="download control light curves in addition to transient light curve",
    )
    parser.add_argument(
        "-n",
        "--num_controls",
        default=None,
        type=int,
        help="number of control light curves per SN",
    )
    parser.add_argument(
        "-r",
        "--radius",
        default=None,
        type=float,
        help="radius of control light curve circle pattern around SN",
    )
    parser.add_argument(
        "--ctrl_coords_filepath",
        type=str,
        default=None,
        help="file name of text file containing table of control light curve coordinates",
    )

    # for downloading single SN with control light curves only
    parser.add_argument(
        "--center_coords",
        type=parse_arg_coords,
        default=None,
        help="comma-separated RA and Dec coordinates of a nearby bright object interfering with the light curve to become center of control light curve circle",
    )

    return parser


if __name__ == "__main__":
    args = define_args().parse_args()
    config = load_config(args.config_file)

    # colnames = load_preset_column_names_from_config("atlas", config)

    # set up directories
    input_dir = config["dir"]["atclean_input"]
    output_dir = config["dir"]["output"]

    # set up Credentials
    creds = Credentials(
        config["credentials"]["atlas_username"],
        config["credentials"]["atlas_password"],
        config["credentials"]["tns_api_key"],
        config["credentials"]["tns_id"],
        config["credentials"]["tns_bot_name"],
    )
    print(f"\nATLAS username: {creds.atlas_username}")
    creds.prompt_for_atlas_password()
    print(f"TNS ID: {creds.tns_id}")
    print(f"TNS bot name: {creds.tns_bot_name}")

    # set up SnInfoTable
    print()
    sninfo_filename = args.sninfo_file or config["dir"]["sninfo_filename"]
    sninfo = SnInfoTable(output_dir, filename=sninfo_filename)

    if len(args.tnsnames) > 1 and (args.sn_coords is not None or args.mjd0 is not None):
        raise ValueError("Cannot apply coordinates and MJD0 to more than one TNS name")

    # parse control light curve args
    num_controls = (
        args.num_controls
        if args.num_controls is not None
        else int(config["download"]["num_controls"])
    )
    radius = (
        args.radius if args.radius is not None else float(config["download"]["radius"])
    )
    sn_min_dist = float(config["download"]["sn_min_dist"])

    downloader = AtlasLightCurveDownloader(
        input_dir, creds.atlas_username, creds.atlas_password
    )

    for tnsname in args.tnsnames:
        print(f"\n--- Downloading ATLAS light curves for {tnsname} ---")

        make_dir_if_not_exists(os.path.join(input_dir, tnsname))
        make_dir_if_not_exists(os.path.join(output_dir, tnsname))

        # get SN location RA and Dec, center location RA and Dec, and MJD0 from command line, SnInfoTable, or TNS API
        sn_coords, center_coords, mjd0 = resolve_sn_coords_and_mjd0(
            tnsname,
            sninfo=sninfo,
            creds=creds,
            arg_mjd0=args.mjd0,
            arg_sn_coords=args.sn_coords,
            arg_center_coords=args.center_coords,
        )
        # update SnInfoTable with resolved RA, Dec, and MJD0
        sninfo.update_row(
            tnsname, coords=sn_coords, mjd0=mjd0, overwrite=args.overwrite
        )
        sninfo.save()

        # construct new ControlCoordinatesTable
        control_coords_table = ControlCoordinatesTableFactory.new(
            args.download_controls,
            num_controls=num_controls,
            radius=radius,
            center_coords=center_coords,
            sn_min_dist=sn_min_dist,
        )

        # download SN and control light curves
        control_coords_table = downloader.download_and_save(
            tnsname,
            sn_coords,
            control_coords_table,
            lookbacktime=args.lookbacktime,
            max_mjd=args.max_mjd,
            overwrite=args.overwrite,
        )
