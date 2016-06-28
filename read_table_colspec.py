#!/usr/bin/env python3
#
# Aaron LI
# Created: 2016-06-25
# Updated: 2016-06-26
#
# Change logs:
# 2016-06-26:
#   * Allow missing data columns
#

"""
Read a table from a plain text file according to a data column
specification.

The "column" here means the column number by each character.
While the "data column" means a meaningful data which stretch multiple
columns. e.g,
------------------------------------------------------------
1234567890123456789012345678901234567890  -> column
  x  xxx.xx  xxx.x   xxxx
 xx   xx.x    xx.xx  xxxx
  |     |       |      |
 DC1   DC2     DC3    DC4                 -> data column
------------------------------------------------------------

This table syntax is usually used in optical astronomy, e.g, SDSS tables.


Sample configuration file:
------------------------------------------------------------
## Configuration for `read_table_colspec.py`
## Date: 2016-06-25

# input table file
infile = <INPUT.TXT>

# output file in CSV format
outfile = <OUTPUT.CSV>

# data column specification
# column number is 1-based, and end column is inclusive.
[colspec]
    # name = col_begin, col_end, type, comment
    ID = 1, 4, int
    RA = 5, 15, float, deg
    Dec = 16, 26, float, deg
    redshift = 27, 34, float
    name = 35, -1, str
------------------------------------------------------------
"""

import re
import argparse
import csv
from pydoc import locate

from configobj import ConfigObj


def parse_colspec(config):
    """
    Parse the data column specification from the config file.

    Parsed colspec syntax:
    [
        (name, col_begin, col_end, type, comment),
        ...
    ]

    NOTE: the above `col_begin` and `col_end` are convert to be 0-based.
    """
    colspec = []
    for name, spec in config.items():
        col_begin, col_end = int(spec[0]), int(spec[1])
        # Convert column number to be 0-based
        col_begin -= 1
        if col_end != -1:
            col_end -= 1
        # Cast from string to type
        # Credit: https://stackoverflow.com/a/29831586/4856091
        t = locate(spec[2])
        if len(spec) == 4:
            comment = spec[3]
        else:
            comment = ""
        colspec.append((name, col_begin, col_end, t, comment))
    return colspec


def parse_line(line, colspec):
    """
    Parse the given line according to the data column specification.

    Parse value syntax:
    [ (name, value, comment), ... ]
    """
    items = []
    for name, col_begin, col_end, t, comment in colspec:
        if col_end == -1:
            value = line[col_begin:].strip()
        else:
            value = line[col_begin:(col_end+1)].strip()
        try:
            value = t(value)
        except ValueError:
            value = None
        items.append((name, value, comment))
    return items


def main():
    parser = argparse.ArgumentParser(
            description="Read table data by column specification")
    parser.add_argument("config", nargs=1,
                        help="config of input, output, column specification")
    args = parser.parse_args()

    config = ConfigObj(args.config[0])
    colspec = parse_colspec(config["colspec"])
    # output column header
    header = [spec[0] if spec[-1] == "" else "%s[%s]" % (spec[0], spec[-1])
              for spec in colspec]
    with open(config["outfile"], "w") as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(header)
        # data lines
        for line in open(config["infile"]).read().splitlines():
            if re.match(r"^\s*$", line):
                # ignore blank lines
                continue
            if re.match(r"^\s*#.*$", line):
                continue
            #
            items = parse_line(line, colspec)
            values = [x[1] for x in items]
            csv_writer.writerow(values)


if __name__ == "__main__":
    main()
