#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Prepare the configuration file for the `sbp_fit.py`.
# And extract name, obsid, r500 information from the '*_INFO.json' file
# to fill the config.
#
# Aaron LI
# 2016-04-21
#
# ChangeLog:
#

import sys
import glob
import os
import re
import json
import argparse
from datetime import datetime


def update_sbpfit_conf(sbpfit_conf, info):
    """
    Update the sbpfit configuration according to the INFO.

    Arguments:
      * sbpfit_conf: list of lines of the sample sbpfit config
      * info: INFO dictionary

    Return:
      updated `sbpfit_conf`
    """
    name = info["Source Name"]
    obsid = int(info["Obs. ID"])

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

    sbpfit_conf_new = []
    for line in sbpfit_conf:
        line_new = re.sub(r"<DATE>", datetime.utcnow().isoformat(), line)
        line_new = re.sub(r"<NAME>", name, line_new)
        line_new = re.sub(r"<OBSID>", "%s" % obsid, line_new)
        line_new = re.sub(r"<R500_PIX>", "%.2f" % r500_pix, line_new)
        line_new = re.sub(r"<R500_KPC>", "%.2f" % r500_kpc, line_new)
        sbpfit_conf_new.append(line_new)

    return sbpfit_conf_new


def main():
    parser = argparse.ArgumentParser(description="Prepare sbpfit config")
    parser.add_argument("-j", "--json", dest="json", required=False,
            help="the *_INFO.json file (default: find ../*_INFO.json)")
    parser.add_argument("-c", "--config", dest="config", required=True,
            help="sample sbpfit configuration")
    parser.add_argument("outfile", nargs="?",
            help="filename of the output sbpfit config " + \
                 "(default: same as the sample config)")
    args = parser.parse_args()

    # default "*_INFO.json"
    info_json = glob.glob("../*_INFO.json")[0]
    if args.json:
        info_json = args.json

    json_str = open(info_json).read().rstrip().rstrip(",")
    info = json.loads(json_str)

    # sample config file
    sbpfit_conf = open(args.config).readlines()

    # output config file
    if args.outfile:
        outfile = args.outfile
    else:
        outfile = os.path.basename(args.config)

    sbpfit_conf_new = update_sbpfit_conf(sbpfit_conf, info)
    with open(outfile, "w") as outf:
        outf.write("".join(sbpfit_conf_new))


if __name__ == "__main__":
    main()

