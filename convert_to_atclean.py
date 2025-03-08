#!/usr/bin/env python

"""
Convert existing non-ATLAS files into ATClean-readable files.

Worry about:
- column names
- possible filters
- filename formatting
"""

from download import load_config
from lightcurve import LightCurve

if __name__ == "__main__":
    args = define_args().parse_args()
    config = load_config(args.config_file)
