#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Extract the excess value and ratio results.
#
# Aaron LI
# Created: 2016-04-27
#

import sys
import os
import json
import csv
import argparse
from collections import OrderedDict


def extract_excess(data):
    """
    Extract the excess value and ratio as well as other misc information
    from the "data" dictionary.
    """
    results = OrderedDict([
            ("name",             data["name"]),
            ("obsid",            data["obsid"]),
            ("model",            data["model"]),
            ("brightness_obs",   data["brightness_obs"]),
            ("brightness_model", data["brightness_model"]),
            ("excess_value",     data["excess_value"]),
            ("excess_ratio",     data["excess_ratio"]),
    ])
    return results


def main():
    parser = argparse.ArgumentParser(
            description="Extract excess results from excess.json")
    parser.add_argument("excess_json", nargs="?", default="excess.json",
            help="json file conatins excess results (default: excess.json)")
    args = parser.parse_args()

    excess_data = json.load(open(args.excess_json))
    results = extract_excess(excess_data)

    # output results
    csv_writer = csv.writer(sys.stdout)
    csv_writer.writerow(results.keys())
    csv_writer.writerow(results.values())


if __name__ == "__main__":
    main()

