#!/usr/bin/env python3
#
# Make the configuration file for `calc_coolfunc.py`
#
# Aaron LI
# Created: 2016-04-21
# Updated: 2016-07-12
#

import glob
import argparse
import subprocess

from info import get_redshift, get_nh


sample_config = """
## Configuration file for `calc_coolfunc.py`

# redshift of the object
redshift = %(redshift).6f

# H column density (unit: 10^22 cm^-2)
nh = %(nh)f

# average abundance (unit: solar)
abundance = %(abund)f
"""


def make_config(info, abund):
    cfg_data = {
        "abund": abund,
        "nh": get_nh(info),
        "redshift": get_redshift(info),
    }
    cfg = sample_config % cfg_data
    return cfg


def get_abund(infile):
    """
    Extract the average abundance from the input file.
    """
    results = subprocess.run(["grep", "^abund", infile],
                             check=True, stdout=subprocess.PIPE)
    abund = float(results.stdout.decode("utf-8").split()[1])
    return abund


def main():
    parser = argparse.ArgumentParser(
        description="Make the config for 'fit_sbp.py'")
    parser.add_argument("-j", "--json", dest="json", required=False,
                        help="the *_INFO.json file " +
                             "(default: find ../*_INFO.json)")
    parser.add_argument("-c", "--cfg", dest="cfg",
                        required=False, default="global.cfg",
                        help="previous config file contains 'abund' " +
                             "(default: global.cfg)")
    parser.add_argument("outfile", nargs="?", default="coolfunc.conf",
                        help="filename of the output coolfunc config " +
                             "(default: coolfunc.conf)")
    args = parser.parse_args()

    # default "*_INFO.json"
    info_json = glob.glob("../*_INFO.json")[0]
    if args.json:
        info_json = args.json

    abund = get_abund(args.cfg)
    config = make_config(info_json, abund)
    open(args.outfile, "w").write(config)


if __name__ == "__main__":
    main()
