#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Extract R500 from the '*_INFO.json' file, and center coordinate from
# existing "sbprofile.reg", and then make the circle regions w.r.t R500
# in order to visualize the FoV coverage of the observations.
#
# Aaron LI
# Created: 2016-04-15
# Updated: 2016-05-04
#
# Changelog:
# 2016-05-04:
#   * Fix a wrong variable
#   * PEP8 fixes
# 2016-04-28:
#   * Add functions "get_r500()" and "get_center()" for code reuse
# 2016-04-16:
#   * Fix R500 unit (convert from kpc to Chandra ACIS pixel)
#

import glob
import re
import json
import argparse


def get_r500(info):
    """
    Get the R500 value (in unit pixel and kpc), as well as the value of
    "kpc_per_pix"
    """
    if isinstance(info, str):
        json_str = open(info).read().rstrip().rstrip(",")
        info = json.loads(json_str)

    if "R500 (kpc)" in info.keys():
        # LWT's
        r500_kpc = float(info["R500 (kpc)"])
    elif "R500" in info.keys():
        # ZZH's
        r500_kpc = float(info["R500"])
    else:
        raise ValueError("Cannot get R500")

    # Convert kpc to Chandra ACIS pixel
    rmax_sbp_pix = float(info["Rmax_SBP (pixel)"])
    rmax_sbp_kpc = float(info["Rmax_SBP (kpc)"])
    kpc_per_pix = rmax_sbp_kpc / rmax_sbp_pix
    r500_pix = r500_kpc / kpc_per_pix

    results = {
            "r500_kpc":    r500_kpc,
            "r500_pix":    r500_pix,
            "kpc_per_pix": kpc_per_pix,
    }
    return results


def get_center(regfile):
    """
    Get the center coordinate from the given region file
    """
    regions = open(regfile).readlines()
    regions = list(filter(lambda x: re.match(r"^(circle|annulus|pie).*",
                                             x, re.I),
                          regions))
    xc, yc = map(float, re.split(r"[(,)]", regions[0])[1:3])
    return (xc, yc)


def frange(x, y, step):
    while x < y:
        yield x
        x += step


def main():
    parser = argparse.ArgumentParser(description="Make R500 circle regions")
    parser.add_argument("-j", "--json", dest="json", required=False,
                        help="the *_INFO.json file " +
                             "(default: find ../*_INFO.json)")
    parser.add_argument("-i", "--regin", dest="regin",
                        required=False, default="sbprofile.reg",
                        help="region from which to extract the center " +
                             "coordinate (default: sbprofile.reg)")
    parser.add_argument("-o", "--regout", dest="regout",
                        required=False, default="r500.reg",
                        help="output region filename (default: r500.reg)")
    args = parser.parse_args()

    # default "*_INFO.json"
    info_json = glob.glob("../*_INFO.json")[0]
    if args.json:
        info_json = args.json

    r500 = get_r500(info_json)
    r500_kpc = r500["r500_kpc"]
    r500_pix = r500["r500_pix"]
    print("R500: %.2f (kpc), %.2f (pixel)" % (r500_kpc, r500_pix))

    # get center coordinate
    xc, yc = get_center(args.regin)

    # output region
    r500_regions = [
            '# Region file format: DS9 version 4.1',
            'global color=white width=1 font="helvetica 10 normal roman"',
            'physical',
    ]
    region_fmt = "circle(%s,%s,{0:.7f}) " % (xc, yc) + \
                 "# text={{{1:.1f} R500 = {2:.1f} pix}}"
    r500_regions.extend([region_fmt.format(r500_pix*ratio, ratio,
                                           r500_pix*ratio)
                         for ratio in frange(0.1, 1.0, 0.1)])
    with open(args.regout, "w") as outfile:
        outfile.write("\n".join(r500_regions) + "\n")


if __name__ == "__main__":
    main()
