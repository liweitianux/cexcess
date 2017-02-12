#!/usr/bin/env python3
#
# Copyright (c) 2016 Aaron LI
# MIT license
#
# Make the configuration file for `deproject_sbp.py`
#
# Created: 2016-07-13
#

import glob
import argparse

from info import get_redshift


sample_config = """
## Configuration file for `deproject_sbp.py`

# redshift of the object (for pixel to distance conversion)
redshift = %(redshift).6f

## SBP extrapolation
# cut radius from which the SBP is fitted for extrapolation,
# specified by the ratio w.r.t sbpfit rc (default: 1.2 * rc)
sbpexp_rignore_ratio = 1.2
# or directly specify the ignorance radius (override above) (unit: pixel)
#sbpexp_rignore = <RIGNORE>
# cut radius to which stop the extrapolation (unit: kpc)
sbpexp_rcut = 3000
"""


def make_config(info):
    cfg_data = {
        "redshift": get_redshift(info),
    }
    cfg = sample_config % cfg_data
    return cfg


def main():
    parser = argparse.ArgumentParser(
        description="Make the config for 'deproject_sbp.py'")
    parser.add_argument("-j", "--json", dest="json", required=False,
                        help="the *_INFO.json file " +
                             "(default: find ../*_INFO.json)")
    parser.add_argument("outfile", nargs="?", default="sbpdeproj.conf",
                        help="filename of the output config " +
                             "(default: sbpdeproj.conf)")
    args = parser.parse_args()

    # default "*_INFO.json"
    info_json = glob.glob("../*_INFO.json")[0]
    if args.json:
        info_json = args.json

    config = make_config(info_json)
    open(args.outfile, "w").write(config)


if __name__ == "__main__":
    main()
