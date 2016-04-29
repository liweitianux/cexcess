#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Aaron LI
# Created: 2016-04-29
# Updated: 2016-04-29
#
# TODO:
#   * to calculate the PEI error
#   * to save a plot with the rectangular and interpolation marked
#

"""
Calculate the power excess index (PEI), which is defined the area ratio of
the lower-left part with respect to the total rectangulr, which is further
defined by the radial power spectrum and the scale of 0.035R500 and 0.35R500,
in the logarithmic space.

Reference:
Zhang, C., et al. 2016, ApJ
"""

__version__ = "0.1.0"
__date__    = "2016-04-29"


import sys
import glob
import argparse
import json
from collections import OrderedDict

import numpy as np
import scipy as sp
import scipy.interpolate
import scipy.integrate

from make_r500_regions import get_r500


def calc_pei(data, r500, interp_np=101):
    """
    Calculate the power excess index (PEI), which is defined the area ratio
    of the lower-left part with respect to the total rectangulr.

    Arguments:
      * data: 3-column power spectrum data (frequency, power, power_err)
      * r500: R500 value in unit of the inverse of the above "frequency"
      * interp_np: number of data points interpolated to calculate PEI
    """
    freqs     = data[:, 0]
    psd1d     = data[:, 1]
    psd1d_err = data[:, 2]
    # frequency values corresponding to 0.35R500 and 0.035R500
    f1 = 1.0 / (0.350 * r500)
    f2 = 1.0 / (0.035 * r500)
    # switch to the logarithmic scale
    # XXX: how to deal with the errors
    mask = (freqs > 0.0)
    x  = np.log10(freqs[mask])
    y  = np.log10(psd1d[mask])
    x1 = np.log10(f1)
    x2 = np.log10(f2)
    # interpolate the power spectrum
    f_interp = sp.interpolate.interp1d(x, y, kind="cubic", assume_sorted=True)
    y1       = f_interp(x1)
    y2       = f_interp(x2)
    if interp_np % 2 == 0:
        # Simpson's rule requires an even number of intervals
        interp_np += 1
    x_interp = np.linspace(x1, x2, num=interp_np)
    y_interp = f_interp(x_interp)
    # calculate the PEI
    area_total = abs(x1 - x2) * abs(y1 - y2)
    area_below = sp.integrate.simps((y_interp-y2), x_interp)
    pei_value  = area_below / area_total
    results = {
            "area_total": area_total,
            "area_below": area_below,
            "pei_value":  pei_value,
            "pei_err":    None,
    }
    return results


def main():
    # default arguments
    default_infile   = "psd.txt"
    default_outfile  = "pei.json"
    default_infojson = "../*_INFO.json"

    parser = argparse.ArgumentParser(
            description="Calculate the power excess index (PEI)",
            epilog="Version: %s (%s)" % (__version__, __date__))
    parser.add_argument("-V", "--version", action="version",
            version="%(prog)s " + "%s (%s)" % (__version__, __date__))
    parser.add_argument("-j", "--json", dest="json", required=False,
            help="the *_INFO.json file (default: find %s)" % default_infojson)
    parser.add_argument("-i", "--infile", dest="infile",
            required=False, default=default_infile,
            help="input data of the radial power spectrum " + \
                 "(default: %s)" % default_infile)
    parser.add_argument("-o", "--outfile", dest="outfile",
            required=False, default=default_outfile,
            help="output json file to save the results " + \
                 "(default: %s)" % default_outfile)
    args = parser.parse_args()

    info_json = glob.glob(default_infojson)[0]
    if args.json:
        info_json = args.json

    json_str    = open(info_json).read().rstrip().rstrip(",")
    info        = json.loads(json_str)
    name        = info["Source Name"]
    obsid       = int(info["Obs. ID"])
    r500        = get_r500(info)
    r500_kpc    = r500["r500_kpc"]
    r500_pix    = r500["r500_pix"]
    kpc_per_pix = r500["kpc_per_pix"]

    psd_data = np.loadtxt(args.infile)
    pei = calc_pei(psd_data, r500=r500_pix)

    pei_data = OrderedDict([
            ("name",            name),
            ("obsid",           obsid),
            ("r500_kpc",        r500_kpc),
            ("r500_pix",        r500_pix),
            ("kpc_per_pix",     kpc_per_pix),
            ("area_total",      pei["area_total"]),
            ("area_below",      pei["area_below"]),
            ("pei",             pei["pei_value"]),
            ("pei_err",         pei["pei_err"]),
    ])
    pei_data_json = json.dumps(pei_data, indent=2)
    print(pei_data_json)
    open(args.outfile, "w").write(pei_data_json+"\n")


if __name__ == "__main__":
    main()

#  vim: set ts=4 sw=4 tw=0 fenc=utf-8 ft=python: #
