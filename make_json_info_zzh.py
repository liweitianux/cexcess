#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Extract the source information and save as a JSON file for each
# source of ZZH from the CSV collection file.
#
# Aaron LI
# 2016-04-14
#

import os
import re
import sys
import csv
import json
import functools
from collections import OrderedDict

argc = len(sys.argv)
if argc < 2:
    print("usage:")
    print("    %s <input_csv> [ <output_json> <obs_id> ]" % \
            os.path.basename(sys.argv[0]))
    sys.exit(1)

# extract default source name and obsid from the output path
cwd = os.getcwd()
m = re.match(r".*zzh/(?P<name>[^_]+)_oi(?P<obsid>\d+)/repro.*$", cwd)
name = m.group("name")
obsid = m.group("obsid")
json_fn = name + "_INFO.json"

csv_fn = sys.argv[1]
if argc >= 3:
    json_fn = sys.argv[2]
if argc == 4:
    obsid = int(sys.argv[3])

with open(csv_fn, "r") as csvfile:
    csv_reader = csv.reader(csvfile)
    csv_data = list(csv_reader)

csv_header = csv_data[0]
obsid_colidx = functools.reduce(None,
        filter(lambda x: x[1] == "Obs. ID", enumerate(csv_header)))[0]
csv_obsid = [ x[obsid_colidx] for x in csv_data ]
obsid_dict = { obsid:idx for idx, obsid in enumerate(csv_obsid) }

obsid_data = csv_data[obsid_dict["%s" % obsid]]
json_data = OrderedDict(zip(csv_header, obsid_data))
with open(json_fn, "w") as jsonfile:
    jsonfile.write(json.dumps(json_data, indent=4, ensure_ascii=False))

#  vim: set ts=4 sw=4 tw=0 fenc=utf-8 ft=python: #
