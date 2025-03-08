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
import pandas as pd
from configparser import ConfigParser
from typing import Dict, List
from download import load_config
from lightcurve import LightCurve, Supernova, SnInfoTable, get_filename
from pdastro import pdastrostatsclass


def parse_config_value(value: str | None):
    # convert 'None' string to actual None
    if value == "None":
        return None
    return value


class PresetColumnNames:
    def __init__(self, config: ConfigParser, preset: str):
        if preset not in config["convert"]:
            raise RuntimeError(f"ERROR: Preset '{preset}' not found in config file.")

        self.preset = preset
        config_preset_settings: Dict[str, str] = config["convert"][self.preset]

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


class ConvertLoop:
    def __init__(
        self,
        preset_colnames: PresetColumnNames,
        filenames: List[str],
        control_indices: List[int],
        input_dir: str,
        output_dir: str,
        sninfo_filename: str = None,
    ):
        self.preset_colnames: PresetColumnNames = preset_colnames

        if len(filenames) < 1:
            raise RuntimeError(
                "ERROR: Please provide at least one file name using the -f argument"
            )
        if len(filenames) != len(control_indices):
            raise RuntimeError(
                f"ERROR: Each file name must have a corresponding control index \n\tfile names (len {len(args.filenames)}): {args.filenames}\n\tcontrol indices (len {len(args.control_indices)}): {args.control_indices}"
            )
        self.filenames: List[str] = filenames

        for control_index in control_indices:
            if not isinstance(control_index, int):  # check if the value is an integer
                raise RuntimeError(
                    f"ERROR: Control index '{control_index}' is not an integer"
                )
            if control_index < 0:  # check if the integer is negative
                raise RuntimeError(
                    f"ERROR: Control index '{control_index}' cannot be negative"
                )
        self.control_indices: List[int] = control_indices

        self.input_dir: str = input_dir
        self.output_dir: str = output_dir
        self.sninfo: SnInfoTable = SnInfoTable(
            self.output_dir, filename=sninfo_filename
        )

    def loop(self, obj_name, ra, dec, mjd0):
        # create ControlCoordinatesTable if needed

        for i in range(len(self.filenames)):
            old_filename = self.filenames[i]
            control_index = self.control_indices[i]

            # load table (.txt or .csv)
            old_lc = pdastrostatsclass()
            if old_filename.endswith(".csv"):
                old_lc.t = pd.read_csv(old_filename)
            else:
                old_lc.load_spacesep(old_filename)

            # move MJD, flux, and dflux columns to the front of the file

            # get filters from filter column
            # if no filter column, set filter to preset
            filts: List[str] = [self.preset_colnames.preset]
            if not self.preset_colnames.filt is None:
                filts = old_lc.t[self.preset_colnames.filt].unique().tolist()

            # handle ra and dec from args and columns

            # add new row to ControlCoordinatesTable if control_index != 0

            # add new row to SnInfoTable if control_index == 0

            # for each filter, save a separate light curve
            for filt in filts:
                new_filename = get_filename(
                    self.output_dir, obj_name, filt, control_index
                )

                # save file

        # save ControlCoordinatesTable


# define command line arguments
def define_args(parser=None, usage=None, conflict_handler="resolve"):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

    parser.add_argument("obj_name", help="transient name")
    parser.add_argument(
        "-p",
        "--preset",
        type=str,
        default=None,
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


if __name__ == "__main__":
    args = define_args().parse_args()
    config = load_config(args.config_file)

    if args.preset is None:
        raise RuntimeError(
            "ERROR: Please specify the preset name to load from the config file (ex. atlas, rubin, tess)"
        )
    preset_colnames = PresetColumnNames(config, args.preset)
    print(preset_colnames.__str__())

    input_dir = config["dir"]["atclean_input"]
    output_dir = config["dir"]["output"]

    convert = ConvertLoop(
        preset_colnames, args.filenames, args.control_indices, input_dir, output_dir
    )
    convert.loop(args.obj_name, args.ra, args.dec, args.mjd0)
