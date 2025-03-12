#!/usr/bin/env python

"""
Convert existing non-ATLAS files into ATClean-readable files.
- Get MJD0 from command line (--mjd0 argument)
- Get RA and Dec from command line (--ra and --dec arguments) or from config file (ra_column_name and dec_column_name)
- Save transient name and coordinates to SnInfoTable in output_dir
- If more than one file name provided, try to create ControlCoordinatesTable using ra_column_name and dec_column_name of light curves with control_index > 0
- Move MJD, flux, and dflux columns to the front of the file
- If filter_column, parse filter column and separate into different light curves
- Save as new file(s) with new input_dir, transient name, and control index
"""

import argparse
import sys
import numpy as np
import pandas as pd
from configparser import ConfigParser
from typing import Dict, List
from download import ControlCoordinatesTable, load_config
from lightcurve import (
    AorB,
    Coordinates,
    LightCurve,
    SnInfoTable,
    get_filename,
)


def parse_config_value(value: str | None):
    # convert 'None' string to actual None
    if value == "None":
        return None
    return value


class PresetColumnNames:
    def __init__(self, config: ConfigParser, preset: str):
        self.preset = preset

        try:
            config_preset_settings: Dict[str, str] = config[f"convert.{preset}"]
        except:
            raise RuntimeError(
                f"ERROR: Preset '{preset}' (field '{f'convert.{preset}'}') not found in config file."
            )

        self.parse_config_values(config_preset_settings)

    def parse_config_values(self, config_preset_settings):
        # required columns
        self.mjd: str = config_preset_settings["mjd_column_name"]
        self.flux: str = config_preset_settings["flux_column_name"]
        self.uncertainty: str = config_preset_settings["uncertainty_column_name"]

        # optional columns
        self.chisquare: str | None = parse_config_value(
            config_preset_settings["chisquare_column_name"]
        )
        self.filt: str | None = parse_config_value(
            config_preset_settings["filter_column_name"]
        )
        self.mag: str | None = parse_config_value(
            config_preset_settings["mag_column_name"]
        )
        self.dmag: str | None = parse_config_value(
            config_preset_settings["dmag_column_name"]
        )
        self.ra: str | None = parse_config_value(
            config_preset_settings["ra_column_name"]
        )
        self.dec: str | None = parse_config_value(
            config_preset_settings["dec_column_name"]
        )

        # extra columns to copy
        columns_to_copy = parse_config_value(config_preset_settings["columns_to_copy"])
        self.columns_to_copy: List[str] = (
            []
            if columns_to_copy is None
            else [col.strip() for col in columns_to_copy.split(",")]
        )

    def __str__(self):
        column_names = [
            f"mjd column: {self.mjd}",
            f"flux column: {self.flux}",
            f"uncertainty column: {self.uncertainty}",
            f"chi-square column: {self.chisquare}",
            f"filter column: {self.filt}",
            f"mag column: {self.mag}",
            f"dmag column: {self.dmag}",
            f"ra column: {self.ra}",
            f"dec column: {self.dec}",
            f"columns to copy: {', '.join(self.columns_to_copy)}",
        ]
        return "\n".join(column_names)


class ConvertLightCurve(LightCurve):
    def __init__(
        self,
        obj_name: str,
        preset_colnames: PresetColumnNames,
        control_index: int = 0,
    ):
        LightCurve.__init__(self, control_index)
        self.obj_name: str = obj_name
        self.preset_colnames: PresetColumnNames = preset_colnames

    def load_raw_t(self, filename: str):
        print(
            f"\n# Loading raw light curve (control index {self.control_index}) at {filename}..."
        )
        if filename.endswith(".csv"):
            self.t = pd.read_csv(filename)
        else:
            self.load_spacesep(filename)

    def move_required_cols_to_front(self):
        cols_to_front = [
            self.preset_colnames.mjd,
            self.preset_colnames.flux,
            self.preset_colnames.uncertainty,
        ]
        self.t = self.t[
            cols_to_front + [col for col in self.t.columns if col not in cols_to_front]
        ]

    def check_single_value_column(self, col_name: str):
        if col_name and self.t[col_name].nunique() != 1:
            raise RuntimeError(
                f"ERROR: Different values found in {col_name} column (control index {self.control_index})"
            )

    def find_coords_in_t(self) -> Coordinates:
        coords = Coordinates()
        if len(self.t) > 0:
            # if RA and Dec columns present, get coords from there
            if not self.preset_colnames.ra is None:
                self.check_single_value_column(self.preset_colnames.ra)
                coords.set_RA(self.t.loc[0, self.preset_colnames.ra])
            if not self.preset_colnames.dec is None:
                self.check_single_value_column(self.preset_colnames.dec)
                coords.set_Dec(self.t.loc[0, self.preset_colnames.dec])
        return coords

    def get_coords(self, arg_ra=None, arg_dec=None) -> Coordinates:
        print("\nSearching for coordinates in command line or light curve...")

        # try to get coordinates from lc columns
        coords_from_t = self.find_coords_in_t()
        if not coords_from_t.is_empty():
            print(f"Found coordinates in light curve: {coords_from_t.__str__()}")

        if self.control_index == 0:
            # try to get coordinates from command line
            coords_from_cmd = Coordinates(arg_ra, arg_dec)
            if not coords_from_cmd.is_empty():
                print(
                    f"Using default coordinates from command line instead: {coords_from_cmd.__str__()}"
                )
                return coords_from_cmd

        if coords_from_t.is_empty():
            print(f"No coordinates found")
        return coords_from_t

    def _save_single_df(self, input_dir, overwrite=False):
        filename = get_filename(
            input_dir,
            self.obj_name,
            filt=self.preset_colnames.preset,
            control_index=self.control_index,
        )
        print(
            f"Saving converted light curve (control index {self.control_index}) with filter {self.preset_colnames.preset}..."
        )
        self.save_lc_by_filename(filename, overwrite=overwrite)

    def _save_dfs_by_filter(self, input_dir, filts, overwrite=False):
        for filt in filts:
            filename = get_filename(
                input_dir,
                self.obj_name,
                filt=filt,
                control_index=self.control_index,
            )
            indices = self.ix_equal(colnames=[self.preset_colnames.filt], val=filt)
            print(
                f"Saving converted light curve (control index {self.control_index}) with filter {filt}..."
            )
            self.save_lc_by_filename(filename, indices=indices, overwrite=overwrite)

    # divide the light curve by filter and save into separate files
    def save(self, input_dir, overwrite=False):
        total_len = len(self.t)
        filt_lens = {}
        if (
            self.preset_colnames.filt is None
        ):  # if no filter column, set filter to preset
            filts = [self.preset_colnames.preset]
        else:  # get filters from filter column
            filts = self.t[self.preset_colnames.filt].unique().tolist()
            for filt in filts:
                filt_lens[filt] = len(
                    np.where(self.t[self.preset_colnames.filt] == filt)[0]
                )

        # sort data by mjd
        self.t = self.t.sort_values(by=[self.preset_colnames.mjd], ignore_index=True)

        # remove rows with duJy=0 or uJy=NaN
        dflux_zero_ix = self.ix_equal(
            colnames=[self.preset_colnames.uncertainty], val=0
        )
        flux_nan_ix = self.ix_is_null(colnames=[self.preset_colnames.flux])
        if len(AorB(dflux_zero_ix, flux_nan_ix)) > 0:
            print(
                f"Deleting {len(dflux_zero_ix) + len(flux_nan_ix)} rows with duJy=0 or uJy=NaN..."
            )
            self.t = self.t.drop(AorB(dflux_zero_ix, flux_nan_ix))

        # save
        if self.preset_colnames.filt is None:
            self._save_single_df(input_dir, overwrite=overwrite)
        else:
            # divide df by filter
            self._save_dfs_by_filter(input_dir, filts, overwrite=overwrite)

        return total_len, filt_lens


class ConvertLoop:
    def __init__(
        self,
        preset_colnames: PresetColumnNames,
        input_dir: str,
        output_dir: str,
        sninfo_filename: str = None,
    ):
        self.preset_colnames: PresetColumnNames = preset_colnames
        self.input_dir: str = input_dir
        self.output_dir: str = output_dir

        print()
        self.sninfo: SnInfoTable = SnInfoTable(
            self.output_dir, filename=sninfo_filename
        )

    def check_args(self, filenames, control_indices):
        if len(filenames) < 1:
            raise RuntimeError(
                "ERROR: Please provide at least one file name using the -f argument"
            )
        if len(filenames) != len(control_indices):
            raise RuntimeError(
                f"ERROR: Each file name must have a corresponding control index \n\tfile names (len {len(args.filenames)}): {args.filenames}\n\tcontrol indices (len {len(args.control_indices)}): {args.control_indices}"
            )

        for control_index in control_indices:
            if not isinstance(control_index, int) or control_index < 0:
                raise RuntimeError(
                    f"Invalid control index: {control_index}. It must be a non-negative integer."
                )

    def loop(
        self,
        obj_name: str,
        filenames: List[str],
        control_indices: List[int],
        ra: str | None = None,
        dec: str | None = None,
        mjd0: float | None = None,
        overwrite: bool = False,
    ):
        self.check_args(filenames, control_indices)

        # create ControlCoordinatesTable if any control light curves present
        ctrl_coords = None
        if len(filenames) > 1:
            ctrl_coords = ControlCoordinatesTable()
            ctrl_coords.num_controls = len(filenames)

        for i in range(len(filenames)):
            old_filename = filenames[i]
            control_index = control_indices[i]

            lc = ConvertLightCurve(
                obj_name,
                self.preset_colnames,
                control_index=control_index,
            )
            lc.load_raw_t(old_filename)

            # move MJD, flux, and dflux columns to the front of the file
            lc.move_required_cols_to_front()

            # try to get coordinates from either command line or ra/dec columns
            coords = lc.get_coords(ra, dec)

            # if not control lc, add new row to SnInfoTable
            if control_index == 0:
                self.sninfo.update_row(
                    obj_name, coords=coords, mjd0=mjd0, overwrite=True
                )

            # for each filter, save a separate light curve
            print()
            total_len, filt_lens = lc.save(self.input_dir, overwrite=overwrite)

            if not ctrl_coords is None:
                # add new row to ControlCoordinatesTable
                ctrl_coords.add_row(
                    obj_name,
                    control_index,
                    coords,
                    ra_offset=np.nan,
                    dec_offset=np.nan,
                    radius=np.nan,
                    n_detec=total_len,
                    filt_lens=filt_lens,
                )

        if not ctrl_coords is None:
            # save ControlCoordinatesTable
            print()
            ctrl_coords.save(self.input_dir, tnsname=obj_name)

        # save SnInfoTable
        self.sninfo.save()


# define command line arguments
def define_args(parser=None, usage=None, conflict_handler="resolve"):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

    parser.add_argument("obj_name", help="transient name")
    parser.add_argument(
        "-p",
        "--preset",
        type=str,
        default="rubin",
        help="preset name from config file (ex. atlas, rubin, tess)",
    )
    parser.add_argument(
        "--mjd0", type=float, default=None, help="transient start date in MJD"
    )
    parser.add_argument(
        "--ra", type=float, default=None, help="transient right ascension"
    )
    parser.add_argument("--dec", type=float, default=None, help="transient declination")
    parser.add_argument(
        "-f",
        "--filenames",
        nargs="+",
        type=str,
        help="one or more file names in raw_input directory to convert",
    )
    parser.add_argument(
        "-i",
        "--control_indices",
        nargs="+",
        type=int,
        default=[0],
        help="one or more ordered control indices corresponding to each file name in -f",
    )
    parser.add_argument(
        "--config_file",
        default="config.ini",
        type=str,
        help="file name of .ini file with settings for this class",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        default=False,
        action="store_true",
        help="overwrite existing file with same file name",
    )

    return parser


if __name__ == "__main__":
    args = define_args().parse_args()
    config = load_config(args.config_file)
    print("Success")

    if args.preset is None:
        raise RuntimeError(
            "ERROR: Please specify the preset name to load from the config file (ex. atlas, rubin, tess)"
        )

    print("\nLoading preset column names from config.ini...")
    preset_colnames = PresetColumnNames(config, args.preset)
    print(preset_colnames.__str__())
    print("Success")

    input_dir = config["dir"]["atclean_input"]
    output_dir = config["dir"]["output"]

    convert = ConvertLoop(preset_colnames, input_dir, output_dir)
    convert.loop(
        args.obj_name,
        args.filenames,
        args.control_indices,
        args.ra,
        args.dec,
        args.mjd0,
        overwrite=args.overwrite,
    )
