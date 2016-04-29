#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Calculate the surface brightness concentration (i.e., C_{SB}), which
# is an index/indicator of the cool core, and may be defined as:
#   (1) brightness(<=40kpc) / brightness(<=400kpc)
#   (2) brightness(<=0.048R500) / brightness(<=0.45R500)
#
# References:
# [1] Santos, J. S., Rosati, P., Tozi, P., et al. 2008, A&A, 483, 35
#
# Aaron LI
# Created: 2016-04-28
# Updated: 2016-04-29
#
# Changelog:
# 2016-04-29:
#   * Fix order of executing ds9 and print message
#   * Fix C_SB calculation
#   * Add reference
#   * Fix "cuspiness" to "concentration"
#   * Add "name" and "obsid" to results
# 2016-04-28:
#   * Add "csb_type" to results
#

import sys
import os
import glob
import re
import json
import argparse
import subprocess
from collections import OrderedDict

from astropy.io import fits

from make_r500_regions import get_r500, get_center


def make_csb_region(regfile, center, r1, r2):
    """
    Make the regions for C_SB and save.
    """
    regions = [
            "pie(%.2f,%.2f,0,%.2f,0,360)" % (center[0], center[1], r1),
            "pie(%.2f,%.2f,0,%.2f,0,360)" % (center[0], center[1], r2),
    ]
    open(regfile, "w").write("\n".join(regions) + "\n")


def calc_csb(evt, expmap, regfile, r1, r2):
    """
    Calculate the C_SB
    """
    csbfile = os.path.splitext(regfile)[0] + ".fits"
    cmd = "dmextract infile='%s[bin sky=@%s]' " % (evt, regfile) + \
          "outfile=%s exp=%s opt=generic clobber=yes" % (csbfile, expmap)
    subprocess.call("punlearn dmextract", shell=True)
    subprocess.call(cmd, shell=True)
    # read calculate C_SB data from output FITS
    with fits.open(csbfile) as csb_fits:
        csb_s_val = csb_fits["HISTOGRAM"].data["SUR_BRI"]
        csb_s_err = csb_fits["HISTOGRAM"].data["SUR_BRI_ERR"]
    # calculate C_SB and error
    area_ratio = (r2 / r1) ** 2
    csb = csb_s_val[0] / csb_s_val[1] / area_ratio
    csb_err = csb * ((csb_s_err[0] / csb_s_val[0])**2 + \
                     (csb_s_err[1] / csb_s_val[1])**2) ** 0.5
    results = OrderedDict([
            ("csb_s1",     csb_s_val[0]),
            ("csb_s1_err", csb_s_err[0]),
            ("csb_s2",     csb_s_val[1]),
            ("csb_s2_err", csb_s_err[1]),
            ("csb",        csb),
            ("csb_err",    csb_err),
    ])
    return results


def main():
    parser = argparse.ArgumentParser(
            description="Calculate the surface brightness concentration")
    # exclusive argument group for C_SB definition
    grp_csb = parser.add_mutually_exclusive_group(required=True)
    grp_csb.add_argument("-K", "--kpc", dest="kpc", action="store_true",
            help="C_SB = brightness(<=0.048R500) / brightness(<=0.45R500)")
    grp_csb.add_argument("-R", "--r500", dest="r500", action="store_true",
            help="C_SB = brightness(<=40kpc) / brightness(<=400kpc)")
    #
    parser.add_argument("-A", "--no-ask", dest="no_ask", required=False,
            action="store_true", help="do NOT check region and ask")
    parser.add_argument("-j", "--json", dest="json", required=False,
            help="the *_INFO.json file (default: find ../*_INFO.json)")
    parser.add_argument("-r", "--region", dest="region",
            required=False, default="sbprofile.reg",
            help="region from which to extract the center coordinate " + \
                 "(default: sbprofile.reg)")
    parser.add_argument("-i", "--infile", dest="infile", required=True,
            help="input energy-filtered EVT2 used to calculate the C_SB")
    parser.add_argument("-e", "--expmap", dest="expmap", required=True,
            help="exposure map of the input image")
    parser.add_argument("-o", "--outfile", dest="outfile", required=True,
            help="output json file to store the C_SB results")
    #
    args = parser.parse_args()

    # default "*_INFO.json"
    info_json = glob.glob("../*_INFO.json")[0]
    if args.json:
        info_json = args.json

    json_str = open(info_json).read().rstrip().rstrip(",")
    info = json.loads(json_str)

    name  = info["Source Name"]
    obsid = int(info["Obs. ID"])

    r500        = get_r500(info)
    r500_kpc    = r500["r500_kpc"]
    r500_pix    = r500["r500_pix"]
    kpc_per_pix = r500["kpc_per_pix"]
    print("R500: %.2f (kpc), %.2f (pixel)" % (r500_kpc, r500_pix))
    # get center coordinate
    xc, yc = get_center(args.region)

    if args.r500:
        csb_type = "r500"
        r1 = 0.048 * r500_pix
        r2 = 0.450 * r500_pix
    elif args.kpc:
        csb_type = "kpc"
        r1 =  40.0 / kpc_per_pix
        r2 = 400.0 / kpc_per_pix
    else:
        raise ValueError("Unknown C_SB definition")

    # make regions for C_SB
    regfile = os.path.splitext(args.outfile)[0] + ".reg"
    make_csb_region(regfile, center=(xc, yc), r1=r1, r2=r2)
    # check region with DS9
    if args.no_ask == False:
        print("Check the C_SB regions; overwrite the region file " + \
                "'%s' if modified" % regfile, flush=True, file=sys.stderr)
        cmd = "ds9 %s -cmap he " % args.infile + \
              "-regions format ciao -regions %s" % regfile
        subprocess.call(cmd, shell=True)
        ans = input("C_SB regions exceed CCD (No/yes/modified)? ")
        if ans == "" or ans[0] in "nN":
            csb_region_note = "OK"
        elif ans[0] in "yY":
            csb_region_note = "EXCESS"
        elif ans[0] in "mM":
            csb_region_note = "MODIFIED"
        else:
            csb_region_note = "???"
    else:
        csb_region_note = None

    # calculate the C_SB
    csb = calc_csb(args.infile, expmap=args.expmap,
            regfile=regfile, r1=r1, r2=r2)
    csb_data = OrderedDict([
            ("name",        name),
            ("obsid",       obsid),
            ("csb_type",    csb_type),
            ("csb_r1",      r1),
            ("csb_r2",      r2),
    ])
    csb_data.update(csb)
    csb_data["csb_region"] = csb_region_note
    csb_data_json = json.dumps(csb_data, indent=2)
    print(csb_data_json)
    open(args.outfile, "w").write(csb_data_json + "\n")


if __name__ == "__main__":
    main()

