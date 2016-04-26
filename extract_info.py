#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Extract R500 from the '*_INFO.json' file, and center coordinate from
# existing "sbprofile.reg", and then make the circle regions w.r.t R500
# in order to visualize the FoV coverage of the observations.
#
# Aaron LI
# 2016-04-22
#
# ChangeLog:
#

import sys
import glob
import os
import re
import json
import csv
import argparse
from collections import OrderedDict


def extract_info(info):
    """
    Extract the requested information (i.e., keywords) from the info
    dictionary.
    """
    # all the keys whose info to be extracted (and keep the order)
    keys = OrderedDict([
            ("Name",               "Source Name"),
            ("ObsID",              ("Obs. ID", int)),
            ("Detector",           "Detector"),
            ("Exposure_ks",        ("Exposure (ks)", float)),
            ("Exposure_clean_ks",  ("Clean Exposure (ks)", float)),
            ("redshift",           ("redshift", float)),
            ("nH",                 None),
            ("Rmax_SBP_pix",       ("Rmax_SBP (pixel)", float)),
            ("Rmax_SBP_kpc",       ("Rmax_SBP (kpc)", float)),
            ("kpc_per_pix",        None),
            ("R500_kpc",           None),
            ("R500_pix",           None),
            ("R500EL_kpc",         None),
            ("R500EU_kpc",         None),
            ("M500",               None),
            ("M500EL",             None),
            ("M500EU",             None),
            ("L500",               None),
            ("L500E",              None),
            ("R200_kpc",           None),
            ("R200_pix",           None),
            ("R200EL_kpc",         None),
            ("R200EU_kpc",         None),
            ("M200",               None),
            ("M200EL",             None),
            ("M200EU",             None),
            ("L200",               None),
            ("L200E",              None),
            ("Tavg_0.2-0.5R500",   ("T(0.2-0.5 R500)", float)),
            ("TavgEL_0.2-0.5R500", ("T_err_l(0.2-0.5 R500)", float)),
            ("TavgEU_0.2-0.5R500", ("T_err_u(0.2-0.5 R500)", float)),
    ])
    # keys corresponding to lwt's INFO
    keys_lwt = {
            "nH":                 ("nH (10^22 cm^-2)", float),
            "R500_kpc":           ("R500 (kpc)", float),
            "R500EL_kpc":         ("R500_err_lower (1sigma)", float),
            "R500EU_kpc":         ("R500_err_upper (1sigma)", float),
            "M500":               ("M500 (M_sun)", float),
            "M500EL":             ("M500_err_lower (1sigma)", float),
            "M500EU":             ("M500_err_upper (1sigma)", float),
            "L500":               ("L500 (erg/s)", float),
            "L500E":              ("L500_err (1sigma)", float),
            "R200_kpc":           ("R200 (kpc)", float),
            "R200EL_kpc":         ("R200_err_lower (1sigma)", float),
            "R200EU_kpc":         ("R200_err_upper (1sigma)", float),
            "M200":               ("M200 (M_sun)", float),
            "M200EL":             ("M200_err_lower (1sigma)", float),
            "M200EU":             ("M200_err_upper (1sigma)", float),
            "L200":               ("L200 (erg/s)", float),
            "L200E":              ("L200_err (1sigma)", float),
    }
    # keys corresponding to zzh's INFO
    keys_zzh = {
            "nH":                 ("nh", float),
            "R500_kpc":           ("R500", float),
            "R500EL_kpc":         ("R500_err_lower", float),
            "R500EU_kpc":         ("R500_err_upper", float),
            "M500":               ("M500 (solar mass)", float),
            "M500EL":             ("M500_err_lower (1 sigma)", float),
            "M500EU":             ("M500_err_upper (1 sigma)", float),
            "L500":               ("L500 (erg/s)", float),
            "L500E":              ("L500_err (1 sigma)", float),
            "R200_kpc":           ("R200", float),
            "R200EL_kpc":         ("R200_err_lower", float),
            "R200EU_kpc":         ("R200_err_upper", float),
            "M200":               ("M200", float),
            "M200EL":             ("M200_err_lower", float),
            "M200EU":             ("M200_err_upper", float),
            "L200":               ("L200", float),
            "L200E":              ("L200_err", float),
    }

    if "IN_SAMPLE" in info.keys():
        # zzh's INFO
        keys.update(keys_zzh)
    else:
        # lwt's INFO
        keys.update(keys_lwt)

    data = keys.copy()
    for k, v in keys.items():
        if v is None:
            continue
        elif isinstance(v, tuple):
            t = v[1]
            try:
                data[k] = t(info[ v[0] ])
            except (ValueError, TypeError):
                data[k] = None
        else:
            data[k] = info[v]

    data["kpc_per_pix"] = data["Rmax_SBP_kpc"] / data["Rmax_SBP_pix"]
    data["R500_pix"]    = data["R500_kpc"] / data["kpc_per_pix"]
    data["R200_pix"]    = data["R200_kpc"] / data["kpc_per_pix"]

    return data


def main():
    parser = argparse.ArgumentParser(description="Extract data from INFO.json")
    parser.add_argument("-j", "--json", dest="json", required=False,
            help="the *_INFO.json file (default: find *_INFO.json)")
    parser.add_argument("path", nargs="?", default=".",
            help="path to the directory contains the INFO.json")
    args = parser.parse_args()

    # path to the directory contains INFO.json
    if args.path == ".":
        path = os.getcwd()
    else:
        path = args.path
        os.chdir(path)

    # default "*_INFO.json"
    info_json = glob.glob("*_INFO.json")[0]
    if args.json:
        info_json = args.json
    json_str = open(info_json).read().rstrip().rstrip(",")
    info = json.loads(json_str)

    data = extract_info(info)
    data["PATH"] = path
    #print(data, file=sys.stderr)

    # output results
    csv_writer = csv.writer(sys.stdout)
    csv_writer.writerow(data.keys())
    csv_writer.writerow(data.values())


if __name__ == "__main__":
    main()

