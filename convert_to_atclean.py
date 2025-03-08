#!/usr/bin/env python

"""
Convert existing non-ATLAS files into ATClean-readable files.
"""

import argparse
from download import load_config
from lightcurve import LightCurve


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
        "--files",
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

    if len(args.files) != len(args.control_indices):
        raise RuntimeError(
            f"ERROR: Each file name must have a corresponding control index \n\tfile names (len {len(args.files)}): {args.files}\n\tcontrol indices (len {len(args.control_indices)}): {args.control_indices}"
        )

    filenames = []
    raw_input_dir = config[args.preset]["raw_input"]
    for filename in args.files:
        if filename.startswith(raw_input_dir):
            filenames.append(filename[len(raw_input_dir) :])
        else:
            filenames.append(filename)
