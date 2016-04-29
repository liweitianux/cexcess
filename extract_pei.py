#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Extract the power excess index (PEI) results.
#
# Aaron LI
# Created: 2016-04-29
#

import sys
import os
import json
import csv
import argparse
from collections import OrderedDict


def extract_pei(data):
    """
    Extract the power excess index (PEI) results from the "data" dictionary.
    """
    results = OrderedDict([
            ("name",            data["name"]),
            ("obsid",           data["obsid"]),
            ("r500_kpc",        data["r500_kpc"]),
            ("r500_pix",        data["r500_pix"]),
            ("kpc_per_pix",     data["kpc_per_pix"]),
            ("area_total",      data["area_total"]),
            ("area_below",      data["area_below"]),
            ("pei",             data["pei"]),
            ("pei_err",         data["pei_err"]),
    ])
    return results


def main():
    parser = argparse.ArgumentParser(
            description="Extract power excess index (PEI) results")
    parser.add_argument("pei_json", help="json file conatins PEI results")
    args = parser.parse_args()

    pei_data = json.load(open(args.pei_json))
    results = extract_pei(pei_data)

    # output results
    csv_writer = csv.writer(sys.stdout)
    csv_writer.writerow(results.keys())
    csv_writer.writerow(results.values())


if __name__ == "__main__":
    main()

