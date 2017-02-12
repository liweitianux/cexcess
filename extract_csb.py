#!/usr/bin/env python3
#
# Copyright (c) 2016 Aaron LI
# MIT license
#
# Extract the surface brightness concentration (i.e., C_{SB}) results.
#
# Created: 2016-04-29
#

import sys
import json
import csv
import argparse
from collections import OrderedDict


def extract_csb(data):
    """
    Extract the surface brightness concentration results
    from the "data" dictionary.
    """
    results = OrderedDict([
            ("name",            data["name"]),
            ("obsid",           data["obsid"]),
            ("csb_type",        data["csb_type"]),
            ("csb_r1",          data["csb_r1"]),
            ("csb_r2",          data["csb_r2"]),
            ("csb_s1",          data["csb_s1"]),
            ("csb_s1_err",      data["csb_s1_err"]),
            ("csb_s2",          data["csb_s2"]),
            ("csb_s2_err",      data["csb_s2_err"]),
            ("csb",             data["csb"]),
            ("csb_err",         data["csb_err"]),
            ("csb_region",      data["csb_region"]),
    ])
    return results


def main():
    parser = argparse.ArgumentParser(
            description="Extract surface brightness concentration results")
    parser.add_argument("csb_json", help="json file conatins C_SB results")
    args = parser.parse_args()

    csb_data = json.load(open(args.csb_json))
    results = extract_csb(csb_data)

    # output results
    csv_writer = csv.writer(sys.stdout)
    csv_writer.writerow(results.keys())
    csv_writer.writerow(results.values())


if __name__ == "__main__":
    main()
