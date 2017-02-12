#!/usr/bin/env python3
#
# Copyright (c) 2016 Aaron LI
# MIT license
#
# Created: 2016-07-11
#

"""
Make the configuration file for `fit_tprofile.py`.
"""

import argparse
import glob

from info import get_redshift


sample_config = """
# redshift of the object (for pixel to distance conversion)
redshift = %(redshift).6f

[model_params]
  # name = initial, lower, upper, variable (FIXED/False to fix the parameter)
  A    = %(A_val).1f, %(A_min).1f, %(A_max).1f, %(A_vary)s
  n    = %(n_val).1f, %(n_min).1f, %(n_max).1f, %(n_vary)s
  xi   = %(xi_val).1f, %(xi_min).1f, %(xi_max).1f, %(xi_vary)s
  a2   = %(a2_val).1f, %(a2_min).1f, %(a2_max).1g, %(a2_vary)s
  a3   = %(a3_val).1f, %(a3_min).1f, %(a3_max).1f, %(a3_vary)s
  beta = %(beta_val).1f, %(beta_min).1f, %(beta_max).1f, %(beta_vary)s
  T0   = %(T0_val).1f, %(T0_min).1f, %(T0_max).1f, %(T0_vary)s
"""


def parse_wang2012_param(pfile):
    """
    Parse the original configuration for `fit_wang2012_model`.
    """
    params = {}
    for name, value, minimum, maximum, vary in map(str.split,
                                                   open(pfile).readlines()):
        par = [
            name,
            float(value),
            float(minimum),
            float(maximum),
            {"T": True, "F": False}.get(vary),
        ]
        params[name] = par
    return params


def make_config(params, redshift):
    cfg_data = {"redshift": redshift}
    for name, value, minimum, maximum, vary in params.values():
        cfg_data.update({
            "%s_val" % name: value,
            "%s_min" % name: minimum,
            "%s_max" % name: maximum,
            "%s_vary" % name: vary,
        })
    cfg = sample_config % cfg_data
    return cfg


def main():
    parser = argparse.ArgumentParser(
        description="Make configuration for 'fit_tprofile.py'")
    parser.add_argument("-j", "--json", dest="json", required=False,
                        help="the *_INFO.json file " +
                             "(default: find ../*_INFO.json)")
    parser.add_argument("-i", "--infile", dest="infile",
                        default="wang2012_param.txt",
                        help="original wang2012 model parameter file " +
                             "(default: wang2012_param.txt")
    parser.add_argument("-o", "--outfile", dest="outfile",
                        default="tprofile.conf",
                        help="output filename")
    args = parser.parse_args()

    # default "*_INFO.json"
    info_json = glob.glob("../*_INFO.json")[0]
    if args.json:
        info_json = args.json
    redshift = get_redshift(info_json)

    params = parse_wang2012_param(args.infile)
    config = make_config(params, redshift)
    open(args.outfile, "w").write(config)


if __name__ == "__main__":
    main()
