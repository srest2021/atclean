#!/usr/bin/env python
"""
@author: Sofia Rest

Zip all converted/download/cleaned light curves of one or more SNe into a single or multiple files.
"""

import argparse
import os
import sys
from typing import Dict, List
import zipfile
from download import load_config


ALLOWED_EXTENSIONS = {
    ".txt",
    ".pdf",
    ".jpg",
    ".png",
    ".md",
    ".ipynb",
    ".py",
    ".json",
    ".ini",
}


def define_args(parser=None, usage=None, conflict_handler="resolve"):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
    parser.add_argument(
        "tnsnames", nargs="+", help="TNS names of the objects to download from ATLAS"
    )
    parser.add_argument(
        "--config_file",
        default="config.ini",
        type=str,
        help="file name of .ini file with settings for this class",
    )
    parser.add_argument(
        "-b",
        "--bulk",
        default=False,
        action="store_true",
        help="store multiple SN files within one zip file",
    )
    return parser


def get_in_dirnames(tnsname: str, input_dir: str, output_dir: str) -> List[str]:
    """
    Get directories to zip for a single SN by its TNS name and the ATLAS input and output directories.
    """
    return [f"{input_dir}/{tnsname}", f"{output_dir}/{tnsname}"]


def get_out_filename(tnsname, output_dir) -> str:
    """
    Get the file name of the output zip file.
    """
    return f"{output_dir}/{tnsname}.zip"


def is_file_allowed(filename: str) -> bool:
    return any(filename.endswith(extension) for extension in ALLOWED_EXTENSIONS)


def get_files_from(dirname) -> Dict[str, str]:
    all_files = {}
    directory_name = os.path.basename(dirname)

    for root, _, files in os.walk(dirname):
        for file in files:
            original_path = os.path.join(root, file)
            relative_path = os.path.join(
                directory_name, os.path.relpath(original_path, start=dirname)
            )
            all_files[original_path] = relative_path
    return all_files


def get_allowed_files_from(dirname) -> List[str]:
    return filter(
        lambda x: is_file_allowed(x) and not os.path.isdir(os.path.join(dirname, x)),
        os.listdir(dirname),
    )


def new_zipfile(out_filename) -> zipfile.ZipFile:
    return zipfile.ZipFile(out_filename, mode="a", compression=zipfile.ZIP_DEFLATED)


def zip_directory(zf: zipfile.ZipFile, in_dirnames: List[str]):
    for in_dirname in in_dirnames:
        filestozip = get_files_from(in_dirname)
        for original_path, relative_path in filestozip.items():
            if not is_file_allowed(relative_path):
                print(f"# Skipping {original_path}")
                continue
            zf.write(original_path, relative_path)
    return zf


def zip_sne_in_bulk(
    tnsnames: List[str], input_dir: str, output_dir: str, out_filename: str
):
    zf = new_zipfile(out_filename)
    for i in range(0, len(tnsnames)):
        print(f"\nZipping {tnsnames[i]} into {out_filename}...")
        in_dirnames = get_in_dirnames(tnsnames[i], input_dir, output_dir)
        zf = zip_directory(zf, in_dirnames)
    zf.close()


def zip_single_sn(
    tnsname: List[str], input_dir: str, output_dir: str, out_filename: str
):
    print(f"\nZipping {tnsname} into {out_filename}...")
    in_dirnames = get_in_dirnames(tnsname, input_dir, output_dir)

    zf = new_zipfile(out_filename)
    zf = zip_directory(zf, in_dirnames)
    zf.close()


if __name__ == "__main__":
    args = define_args().parse_args()

    tnsnames = args.tnsnames
    bulk = args.bulk
    print(
        f"List of transients to zip {'in bulk' if bulk else 'individually'}: {tnsnames}"
    )

    config = load_config(args.config_file)

    input_dir = config["dir"]["atclean_input"]
    output_dir = config["dir"]["output"]
    print(f"ATClean input directory: {input_dir}")
    print(f"Output directory: {output_dir}")

    if bulk:
        out_filename = f"{output_dir}/lcs.zip"
        zip_sne_in_bulk(tnsnames, input_dir, output_dir, out_filename)
    else:
        for i in range(len(tnsnames)):
            out_filename = get_out_filename(tnsnames[i], output_dir)
            zip_single_sn(tnsnames[i], input_dir, output_dir, out_filename)
