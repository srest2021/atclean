#!/usr/bin/env python
"""
Unzip light curve files and place them back into their appropriate directories
based on filename pattern.

@author: Sofia Rest
"""

import argparse
import os
import re
import zipfile
from download import load_config


def define_args(parser=None, usage=None, conflict_handler="resolve"):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
    parser.add_argument("filepaths", nargs="+", help="Path(s) to zip file(s) to unzip.")
    parser.add_argument(
        "--config_file",
        default="config.ini",
        type=str,
        help="file name of .ini file with settings for this class",
    )
    return parser


def is_input_file(filename: str) -> bool:
    """
    Determine if a file belongs in the input_dir based on its name.
    Handles:
    - <tnsname>/<tnsname>.<filter>.lc.txt
    - controls/<tnsname>_i<###>.<filter>.lc.txt
    """
    basename = os.path.basename(filename)

    # Match regular input light curves: e.g. 2019vxm.c.lc.txt
    if re.match(r"^[^.]+\.[^.]+\.lc\.txt$", basename):
        return True

    # Match control light curves: e.g. 2019vxm_i002.c.lc.txt
    if re.match(r"^[^.]+_i\d{3}\.[^.]+\.lc\.txt$", basename) and filename.startswith(
        "controls/"
    ):
        return True

    return False


def extract_file(zf: zipfile.ZipFile, member: str, input_dir: str, output_dir: str):
    target_base = input_dir if is_input_file(member) else output_dir
    target_path = os.path.join(target_base, member)

    # Ensure the directory exists
    os.makedirs(os.path.dirname(target_path), exist_ok=True)

    # Extract the file to the correct location
    with zf.open(member) as source, open(target_path, "wb") as target:
        target.write(source.read())

    label = "ATClean input" if target_base == input_dir else "output"
    print(f"✓ {member} → <{label} directory>/{member}")


def unzip_file(path: str, input_dir: str, output_dir: str):
    print(f"\n📦 Unzipping: {os.path.basename(path)}")
    with zipfile.ZipFile(path, "r") as zf:
        for member in zf.namelist():
            extract_file(zf, member, input_dir, output_dir)
    print("✅ Done.\n")


if __name__ == "__main__":
    args = define_args().parse_args()

    config = load_config(args.config_file)
    input_dir = config["dir"]["atclean_input"]
    output_dir = config["dir"]["output"]

    print(f"\n📁 ATClean input directory:  {input_dir}")
    print(f"📁 Output directory: {output_dir}")

    for path in args.filepaths:
        unzip_file(path, input_dir, output_dir)
