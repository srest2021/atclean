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
from typing import List
from download import load_config
from lightcurve import LightCurve, SnInfoTable


class ConvertLoop:
    def __init__(
        self,
        preset_settings,
        filenames: List[str],
        control_indices: List[int],
        input_dir: str,
        output_dir: str,
        sninfo_filename: str = None,
    ):
        self.preset_settings = preset_settings
        self.filenames = filenames
        self.control_indices = control_indices

        self.input_dir = input_dir
        self.output_dir = output_dir
        self.sninfo: SnInfoTable = SnInfoTable(
            self.output_dir, filename=sninfo_filename
        )


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


if __name__ == "__main__":
    args = define_args().parse_args()
    config = load_config(args.config_file)

    if args.preset is None:
        raise RuntimeError(
            "ERROR: Please specify the preset name to load from the config file (ex. atlas, rubin, tess)"
        )
    preset_settings = config["convert"][args.preset]

    if len(args.filenames) < 1:
        raise RuntimeError(
            "ERROR: Please provide at least one file name using the -f argument"
        )
    if len(args.filenames) != len(args.control_indices):
        raise RuntimeError(
            f"ERROR: Each file name must have a corresponding control index \n\tfile names (len {len(args.files)}): {args.files}\n\tcontrol indices (len {len(args.control_indices)}): {args.control_indices}"
        )

    input_dir = config["dir"]["atclean_input"]
    output_dir = config["dir"]["output"]

    convert = ConvertLoop(
        preset_settings, args.filenames, args.control_indices, input_dir, output_dir
    )
