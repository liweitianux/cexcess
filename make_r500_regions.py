#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Extract R500 from the '*_INFO.json' file, and center coordinate from
# existing "sbprofile.reg", and then make the circle regions w.r.t R500
# in order to visualize the FoV coverage of the observations.
#
# Aaron LI
# 2016-04-15
#
# ChangeLog:
# 2016-04-16:
#   * Fix R500 unit (convert from kpc to Chandra ACIS pixel)
#

import sys
import glob
import os
import re
import json
import argparse


def frange(x, y, step):
    while x < y:
        yield x
        x += step


def main():
    parser = argparse.ArgumentParser(description="Make R500 circle regions")
    parser.add_argument("-j", "--json", dest="json", required=False,
            help="the *_INFO.json file (default: find ../*_INFO.json)")
    parser.add_argument("-i", "--regin", dest="regin",
            required=False, default="sbprofile.reg",
            help="region from which to extract the center coordinate " + \
                 "(default: sbprofile.reg)")
    parser.add_argument("-o", "--regout", dest="regout",
            required=False, default="r500.reg",
            help="output region filename (default: r500.reg)")
    args = parser.parse_args()

    # default "*_INFO.json"
    info_json = glob.glob("../*_INFO.json")[0]
    if args.json:
        info_json = args.json

    json_str = open(info_json).read().rstrip().rstrip(",")
    info = json.loads(json_str)

    if "R500 (kpc)" in info.keys():
        # lwt's
        r500_kpc = float(info["R500 (kpc)"])
    elif "R500" in info.keys():
        # zzh's
        r500_kpc = float(info["R500"])
    else:
        raise ValueError("Cannot get R500 from INFO.json")

    # Convert kpc to Chandra ACIS pixel
    rmax_sbp_pix = float(info["Rmax_SBP (pixel)"])
    rmax_sbp_kpc = float(info["Rmax_SBP (kpc)"])
    r500_pix = r500_kpc / rmax_sbp_kpc * rmax_sbp_pix
    print("R500: %.2f (kpc), %.2f (pixel)" % (r500_kpc, r500_pix))

    # get center coordinate
    regions = open(args.regin).readlines()
    regions = list(filter(lambda x: re.match(r"^(circle|annulus|pie).*",
                                             x, re.I),
                          regions))
    xc, yc = re.split(r"[(,)]", regions[0])[1:3]

    # output region
    r500_regions = [
            '# Region file format: DS9 version 4.1',
            'global color=white width=1 font="helvetica 10 normal roman"',
            'physical',
    ]
    region_fmt = "circle(%s,%s,{0:.7f}) # text={{{1:.1f} R500 = {2:.1f} pix}}" \
            % (xc, yc)
    r500_regions.extend([ region_fmt.format(r500_pix*ratio, ratio,
                                            r500_pix*ratio)
                          for ratio in frange(0.1, 1.0, 0.1) ])
    with open(args.regout, "w") as outfile:
        outfile.write("\n".join(r500_regions) + "\n")


if __name__ == "__main__":
    main()

