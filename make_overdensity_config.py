#!/usr/bin/env python3
#
# Copyright (c) 2016 Aaron LI
# MIT license
#
# Make the configuration file for `calc_overdensity.py`
#
# Created: 2016-07-13
#

import glob
import argparse

from info import get_redshift


sample_config = """
## Configuration file for `calc_overdensity.py`

# redshift of the source (critical density)
redshift = %(redshift).6f
"""


def make_config(info):
    cfg_data = {
        "redshift": get_redshift(info),
    }
    cfg = sample_config % cfg_data
    return cfg


def main():
    parser = argparse.ArgumentParser(
        description="Make the config for 'calc_overdensity.py'")
    parser.add_argument("-j", "--json", dest="json", required=False,
                        help="the *_INFO.json file " +
                             "(default: find ../*_INFO.json)")
    parser.add_argument("outfile", nargs="?", default="overdensity.conf",
                        help="filename of the output overdensity config " +
                             "(default: overdensity.conf)")
    args = parser.parse_args()

    # default "*_INFO.json"
    info_json = glob.glob("../*_INFO.json")[0]
    if args.json:
        info_json = args.json

    config = make_config(info_json)
    open(args.outfile, "w").write(config)


if __name__ == "__main__":
    main()
