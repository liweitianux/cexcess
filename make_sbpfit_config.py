#!/usr/bin/env python3
#
# Copyright (c) 2016 Aaron LI
# MIT license
#
# Make the configuration file for the `sbp_fit.py`.
# The name, obsid, r500 information is extracted from the '*_INFO.json' file.
#
# Created: 2016-04-21
#
# Change logs:
# 2016-07-12:
#   * Rewrite according to `make_tprofile_config.py`,
#     and make the sample config built-in
# 2016-04-26:
#   * Minor update to output file write
#

import glob
import argparse

from info import get_name, get_obsid, get_r500


sample_config = """
## Configuration for `fit_sbp.py`

name     = %(name)s
obsid    = %(obsid)d
r500_pix = %(r500_pix)f
r500_kpc = %(r500_kpc)f

# sbp model: "sbeta" or "dbeta"
model    = sbeta

[sbeta]
outfile     = sbpfit_sbeta.json
imgfile     = sbpfit_sbeta.png
#ignore      = 0.0-20.0,
ignore_r500 = 0.0-0.1,
  [[params]]
  # name = initial, lower, upper, variable (FIXED/False to fix the parameter)
  s0    = 1.0e-8,  0.0,  1.0e-6
  rc    = 30.0,    5.0,  1.0e4
  beta  = 0.7,     0.3,  1.1
  bkg   = 1.0e-10, 0.0,  1.0e-8
"""


def make_config(info):
    cfg_data = {
        "name": get_name(info)["name"],
        "obsid": get_obsid(info),
        "r500_pix": get_r500(info)["r500_pix"],
        "r500_kpc": get_r500(info)["r500_kpc"],
    }
    cfg = sample_config % cfg_data
    return cfg


def main():
    parser = argparse.ArgumentParser(
        description="Make the config for 'fit_sbp.py'")
    parser.add_argument("-j", "--json", dest="json", required=False,
                        help="the *_INFO.json file " +
                             "(default: find ../*_INFO.json)")
    parser.add_argument("outfile", nargs="?", default="sbpfit.conf",
                        help="filename of the output sbpfit config " +
                             "(default: sbpfit.conf)")
    args = parser.parse_args()

    # default "*_INFO.json"
    info_json = glob.glob("../*_INFO.json")[0]
    if args.json:
        info_json = args.json

    config = make_config(info_json)
    open(args.outfile, "w").write(config)


if __name__ == "__main__":
    main()
