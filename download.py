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

from typing import Dict, List, Optional, Type
import os, sys, requests, argparse, configparser, math
import pandas as pd
import numpy as np
from getpass import getpass
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from lightcurve import FullLightCurve, LightCurve
from utils import (
    Coordinates,
    Credentials,
    PresetColumnNames,
    SnInfoTable,
    find_all_control_indices,
    find_all_filts,
    is_sn_in_subdir,
    load_config,
    load_preset_column_names_from_config,
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
        self.closebright: bool = False
        self.closebright_min_dist: Optional[float] = None

    def init_load(self, directory: str, tnsname: str):
        filename = self.get_filepath(directory, tnsname)
        print(f"Loading control coordinates table at {filename}...")
        self.load(filename)
        print("Success")
        print(self.__str__())

        self._set_num_controls_from_t()

    def init_read_from_file(self, filepath: str):
        print(f"Loading control coordinates table at {filepath}...")
        self._read(filepath)
        print("Success")
        print(self.__str__())

        self._set_num_controls_from_t()

    def init_closebright(
        self, center_coords: Coordinates, num_controls: int, closebright_min_dist: float
    ):
        self.num_controls = num_controls
        # set self.radius later

        self.closebright = True
        self.center_coords = center_coords
        self.closebright_min_dist = closebright_min_dist

        print(f"Setting circle pattern of {self.num_controls} control light curves")

    def init_default(self, num_controls: int, radius: float):
        self.num_controls = num_controls
        self.radius = Angle(radius, u.arcsec)
        # set self.center_coords later

        print(
            f"Setting circle pattern of {self.num_controls} control light curves with radius of {self.radius} from center"
        )

    def ready_for_construction(self):
        self.reset_table()

        assert self.num_controls is not None
        assert self.radius is not None
        assert self.sn_coords is not None
        assert self.center_coords is not None

        if self.closebright:
            assert self.closebright_min_dist is not None

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

        if self.closebright:
            # check to see if control light curve location is within minimum distance from SN location
            offset_sep = self.sn_coords.get_distance(coords).arcsecond
            if offset_sep < self.closebright_min_dist:
                print(
                    f'Control light curve {control_index:3d} too close to SN location ({offset_sep}" away) with minimum distance to SN as {self.closebright_min_dist}"; skipping control light curve...'
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

    def construct(self, tnsname: str, full_sn_lc: FullLightCurve):
        self.sn_coords = full_sn_lc.coords

        # add row for SN position
        total_len, filt_lens = full_sn_lc.get_filt_lens()
        if self.closebright:
            # circle pattern radius is distance between SN and close bright object
            self.radius = self.sn_coords.get_distance(self.center_coords)

            self.add_row(
                tnsname,
                0,
                self.sn_coords,
                ra_offset=np.nan,
                dec_offset=np.nan,
                radius=self.radius,
                n_detec=total_len,
                filt_lens=filt_lens,
            )
        else:
            # center coordinates are the SN location
            self.center_coords = self.sn_coords

            self.add_row(
                tnsname, 0, self.sn_coords, n_detec=total_len, filt_lens=filt_lens
            )

        self.ready_for_construction()

        # add row for each control light curve
        for i in range(1, self.num_controls + 1):
            self.construct_row(i)

        print("Control light curve coordinates generated: \n", self.__str__())

    def get_filepath(self, directory, tnsname):
        return os.path.join(directory, tnsname, f"{tnsname}_control_coords.txt")

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

        print(f"Saving control coordinates table at {filepath}...")
        if self.t is None:
            raise RuntimeError(
                "Cannot save ControlCoordinatesTable: table (self.t) is None"
            )
        if overwrite or not os.path.exists(filepath):
            self.t.to_string(filepath, index=False)

    def __str__(self):
        if self.t is None:
            return "ControlCoordinatesTable is empty (no data loaded)"
        with pd.option_context("display.float_format", "{:,.8f}".format):
            return self.t.to_string()


class AtlasAuthenticator:
    @staticmethod
    def authenticate(username: str, password: str) -> Dict[str, str]:
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


def parse_arg_coords(arg_coords: str) -> Coordinates:
    parsed_coords = parse_comma_separated_string(arg_coords)
    if parsed_coords is None:
        raise RuntimeError(
            f"Parsing comma-separated --coords argument failed: {arg_coords}"
        )

    if len(parsed_coords) > 2:
        raise RuntimeError(
            "Too many coordinates in --coords argument! Please provide comma-separated RA and Dec onlyy."
        )
    if len(parsed_coords) < 2:
        raise RuntimeError(
            "Too few coordinates in --coords argument! Please provide comma-separated RA and Dec."
        )
    return Coordinates(parsed_coords[0], parsed_coords[1])


class DownloadLoop:
    def __init__(
        self,
        input_dir: str,
        output_dir: str,
        creds: Credentials,
        sninfo: SnInfoTable,
        # controls: ControlCoordinatesTable,
        colnames: PresetColumnNames,
        overwrite: bool = False,
    ):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.overwrite = overwrite

        self.lcs: Dict[int, FullLightCurve] = {}

        self.controls: Optional[ControlCoordinatesTable] = None

        # SN info table
        self.sninfo: SnInfoTable = sninfo

        # ATLAS and TNS credentials
        self.creds: Credentials = creds

        # ATLAS preset column names
        self.colnames: PresetColumnNames = colnames

    def init_control_coordinates_table(
        self,
        ctrl_coords_filepath: Optional[str] = None,
        num_controls: Optional[int] = None,
        radius: Optional[float] = None,
        center_coords: Optional[Coordinates] = None,
        closebright_min_dist: Optional[float] = None,
    ):
        self.controls = ControlCoordinatesTable()

        if ctrl_coords_filepath:
            self.controls.init_read_from_file(ctrl_coords_filepath)

        else:
            if center_coords:
                if num_controls is None:
                    raise ValueError(
                        "`num_controls` must be provided when using closebright strategy."
                    )
                if closebright_min_dist is None:
                    raise ValueError(
                        "`closebright_min_dist` must be provided when using closebright strategy."
                    )

                # TODO: option to parse from SN info table

                self.controls.init_closebright(
                    center_coords,
                    num_controls,
                    closebright_min_dist,
                )
            else:
                if num_controls is None:
                    raise ValueError(
                        "`num_controls` must be provided when using closebright strategy."
                    )
                if radius is None:
                    raise ValueError(
                        "`radius` must be provided when using closebright strategy."
                    )
                self.controls.init_default(num_controls, radius)

    def _construct_full_sn_lc(
        self,
        tnsname,
        arg_coords: Optional[Coordinates] = None,
        arg_mjd0: Optional[float] = None,
    ):
        # first, try SN info table
        ra, dec, mjd0 = self.sninfo.get_info(tnsname)

        # next, overwrite defaults with command line args
        if arg_coords is not None:
            ra, dec = arg_coords.get_RA_str(), arg_coords.get_Dec_str()
            print(f"Setting coordinates to --coords argument: RA {ra}, Dec {dec}")
        if arg_mjd0 is not None:
            mjd0 = arg_mjd0
            print(f"Setting MJD0 to --mjd0 argument: {mjd0} MJD")

        try:
            self.lcs[0] = FullLightCurve(0, ra, dec, mjd0)
        except Exception as e:
            print(
                f"WARNING: Could not construct light curve object with RA {ra}, Dec {dec}, and MJD0 {mjd0}: {str(e)}."
            )
            self.lcs[0] = FullLightCurve(0)

        # try to query TNS for any missing data
        self.creds.validate_tns_credentials()
        self.lcs[0].get_tns_data(
            tnsname,
            self.creds.tns_api_key,
            self.creds.tns_id,
            self.creds.tns_bot_name,
        )

        # add final RA, Dec, MJD0 to SN info table
        self.sninfo.update_row(
            tnsname, self.lcs[0].coords, self.lcs[0].mjd0, overwrite=self.overwrite
        )

    def _load_or_download_sn_lc(
        self,
        headers: Dict[str, str],
        tnsname: str,
        lookbacktime: Optional[float] = None,
        max_mjd: Optional[float] = None,
        download_controls: bool = False,
    ):
        if not self.overwrite and is_sn_in_subdir(self.input_dir, tnsname):
            print(
                f"Overwrite set to False and SN light curve already exists; skipping download..."
            )
            if download_controls and self.controls.t is None:
                self.controls.load(self.input_dir, tnsname)
        else:
            # download SN light curve
            self.lcs[0].download(headers, lookbacktime=lookbacktime, max_mjd=max_mjd)
            self.lcs[0].save(self.colnames, self.input_dir, tnsname, overwrite=True)

            if download_controls and self.controls.t is None:
                self.controls.construct(tnsname, self.lcs[0])

                # save the constructed table without the filt lens
                # so that if it crashes later, we can just load it again
                self.controls.save(self.input_dir, tnsname)

    def _download_control_lightcurves(
        self,
        headers: Dict[str, str],
        tnsname: str,
        lookbacktime: Optional[float] = None,
        max_mjd: Optional[float] = None,
    ):
        existing_control_indices = find_all_control_indices(self.input_dir, tnsname)

        for i in range(1, len(self.controls.t)):
            control_index = self.controls.t.at[i, "control_index"]
            print(f"\nControl light curve {control_index}")

            if not self.overwrite and control_index in existing_control_indices:
                print(
                    f"Overwrite set to {self.overwrite} and light curve already exists; skipping download..."
                )

                # load saved lcs and update filt lens
                filt_lens = {}
                for filt in find_all_filts(self.input_dir, tnsname):
                    lc = LightCurve(
                        self.colnames, control_index=control_index, filt=filt
                    )
                    lc.load_lc(self.input_dir, tnsname)
                    filt_lens[filt] = len(lc.t)
                total_len = sum([filt_len for filt_len in filt_lens.values()])

            else:
                self.lcs[control_index] = FullLightCurve(
                    control_index,
                    self.controls.t.at[i, "ra"],
                    self.controls.t.at[i, "dec"],
                )
                self.lcs[control_index].download(
                    headers, lookbacktime=lookbacktime, max_mjd=max_mjd
                )
                self.lcs[control_index].save(
                    self.colnames, self.input_dir, tnsname, overwrite=self.overwrite
                )

                total_len, filt_lens = self.lcs[control_index].get_filt_lens()

            self.controls.update_filt_lens(control_index, total_len, filt_lens)

        # save control coordinates table
        # always overwrite the previously saved version with no filt lens
        self.controls.save(self.input_dir, tnsname, overwrite=True)

    def download_lcs(
        self,
        headers: Dict[str, str],
        tnsname: str,
        download_controls: bool = False,
        arg_coords: Optional[Coordinates] = None,
        arg_mjd0: Optional[float] = None,
        lookbacktime: Optional[float] = None,
        max_mjd: Optional[float] = None,
    ):
        print(f"\nDOWNLOADING ATLAS LIGHT CURVES FOR: SN {tnsname}\n")
        self.lcs: Dict[int, FullLightCurve] = {}

        try:
            self._construct_full_sn_lc(
                tnsname, arg_coords=arg_coords, arg_mjd0=arg_mjd0
            )
        except Exception as e:
            print(
                f"Could not construct light curve object: {str(e)}. Skipping to next SN..."
            )
            return

        self._load_or_download_sn_lc(
            headers,
            tnsname,
            lookbacktime=lookbacktime,
            max_mjd=max_mjd,
            download_controls=download_controls,
        )

        self.sninfo.save()

        if download_controls:
            self._download_control_lightcurves(
                headers, tnsname, lookbacktime=lookbacktime, max_mjd=max_mjd
            )

    def loop(
        self,
        tnsnames: List[str],
        arg_coords: Optional[Coordinates] = None,
        arg_mjd0: Optional[float] = None,
        lookbacktime: Optional[float] = None,
        max_mjd: Optional[float] = None,
        download_controls: bool = False,
        ctrl_coords_filepath: str | None = None,
        num_controls: int | None = None,
        radius: float | None = None,
        center_coords: Coordinates | None = None,
        closebright_min_dist: float | None = None,
    ):
        if len(tnsnames) > 1 and (arg_coords or arg_mjd0):
            raise ValueError(
                "Cannot apply coordinates and MJD0 to more than one TNS name"
            )

        print("\nConnecting to ATLAS API...")
        headers = AtlasAuthenticator.authenticate(
            self.creds.atlas_username, self.creds.atlas_password
        )

        if headers is None:
            raise RuntimeError("No token header!")

        for obj_index in range(len(tnsnames)):
            if download_controls:
                self.init_control_coordinates_table(
                    ctrl_coords_filepath=ctrl_coords_filepath,
                    num_controls=num_controls,
                    radius=radius,
                    center_coords=center_coords,
                    closebright_min_dist=closebright_min_dist,
                )

            self.download_lcs(
                headers,
                tnsnames[obj_index],
                arg_coords=arg_coords,
                arg_mjd0=arg_mjd0,
                lookbacktime=lookbacktime,
                max_mjd=max_mjd,
                download_controls=download_controls,
            )


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
        type=str,
        default=None,
        help="comma-separated RA and Dec of SN light curve to download",
    )
    parser.add_argument(
        "--mjd0", type=float, default=None, help="transient start date in MJD"
    )

    # for control light curves
    parser.add_argument(
        "-c",
        "--controls",
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
        type=str,
        default=None,
        help="comma-separated RA and Dec coordinates of a nearby bright object interfering with the light curve to become center of control light curve circle",
    )

    return parser


if __name__ == "__main__":
    args = define_args().parse_args()
    config = load_config(args.config_file)

    colnames = load_preset_column_names_from_config("atlas", config)

    # set up directories
    input_dir = config["dir"]["atclean_input"]
    output_dir = config["dir"]["output"]
    make_dir_if_not_exists(input_dir)
    make_dir_if_not_exists(output_dir)

    # set up credentials
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

    if not args.controls and (
        args.ctrl_coords_filepath
        or args.center_coords
        or args.num_controls
        or args.radius
    ):
        raise RuntimeError(
            "Please specify control light curve downloading (-c or --controls) before using any of the following arguments: --ctrl_coords, --closebright, --num_controls, --radius."
        )

    download = DownloadLoop(
        input_dir,
        output_dir,
        creds,
        sninfo,
        colnames,
        overwrite=args.overwrite,
    )

    download.loop(
        args.tnsnames,
        arg_coords=parse_arg_coords(args.sn_coords) if args.sn_coords else None,
        arg_mjd0=args.mjd0,
        lookbacktime=args.lookbacktime,
        max_mjd=args.max_mjd,
        download_controls=args.controls,
        ctrl_coords_filepath=args.ctrl_coords_filepath,
        num_controls=(
            args.num_controls
            if args.num_controls
            else int(config["download"]["num_controls"])
        ),
        radius=args.radius if args.radius else float(config["download"]["radius"]),
        center_coords=args.center_coords,
        closebright_min_dist=float(config["download"]["closebright_min_dist"]),
    )
