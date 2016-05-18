#!/usr/bin/env python3
#
# Read results in JSON format and output as CSV format.
#
# Aaron LI
# Created: 2016-04-27
# Updated: 2016-05-06
#

import sys
import json
import csv
import argparse
from collections import OrderedDict


def main():
    parser = argparse.ArgumentParser(
            description="Read JSON results and output as CSV format")
    parser.add_argument("json", help="input JSON file")
    parser.add_argument("csv", nargs="?", help="optional output CSV file")
    args = parser.parse_args()

    results = json.load(open(args.json), object_pairs_hook=OrderedDict)

    csv_writer = csv.writer(sys.stdout)
    csv_writer.writerow(results.keys())
    csv_writer.writerow(results.values())

    if args.csv:
        with open(args.csv, "w") as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(results.keys())
            csv_writer.writerow(results.values())


if __name__ == "__main__":
    main()
