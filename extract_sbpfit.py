#!/usr/bin/env python3
#
# Copyright (c) 2016 Aaron LI
# MIT license
#
# Extract the sbpfit results, and output in CSV format.
#
# Created: 2016-04-27
#

import sys
import os
import json
import csv
import argparse
from collections import OrderedDict

from configobj import ConfigObj


def extract_sbpfit(data, config):
    """
    Extract the SBP fitting results (sbpfit) as well as some config
    information, and return the results as a Python dictionary.
    """
    results = OrderedDict([
            # basic information from config
            ("name",        config["name"]),
            ("obsid",       int(config["obsid"])),
            ("r500_pix",    float(config["r500_pix"])),
            ("r500_kpc",    float(config["r500_kpc"])),
            ("unit",        config["unit"]),
            ("model",       config["model"]),
            # basic fitting results
            ("ndata",       data["ndata"]),
            ("nvarys",      data["nvarys"]),
            ("nfree",       data["nfree"]),
            ("nfev",        data["nfev"]),
            ("chisqr",      data["chisqr"]),
            ("redchi",      data["redchi"]),
            ("aic",         data["aic"]),
            ("bic",         data["bic"]),
    ])
    # fitted paramters value and confidence intervals
    results_params = extract_params(data)
    results.update(results_params)
    return results


def extract_params(data):
    """
    Extract the values and confidence intervals (if present)
    for each parameter.
    """
    results = OrderedDict()
    for pname, pvalue in data["params"].items():
        # best value
        results[pname] = pvalue[0]
        # confidence intervals (if present)
        if "ci" in data:
            for ci_name, ci_value in data["ci"][pname].items():
                if ci_name == "best":
                    continue
                results["%s_%s_L" % (pname, ci_name)] = ci_value[0]
                results["%s_%s_U" % (pname, ci_name)] = ci_value[1]
    return results


def main():
    parser = argparse.ArgumentParser(
            description="Extract SBP fitting results")
    parser.add_argument("config", nargs="?", default="sbpfit.conf",
            help="config used for SBP fitting (default: sbpfit.conf)")
    args = parser.parse_args()

    config = ConfigObj(args.config)
    modelname = config["model"]
    sbpfit_outfile = config.get("outfile")
    sbpfit_outfile = config[modelname].get("outfile", sbpfit_outfile)
    sbpfit_outfile = os.path.join(os.path.dirname(args.config), sbpfit_outfile)
    sbpfit_results = json.load(open(sbpfit_outfile),
                               object_pairs_hook=OrderedDict)
    results = extract_sbpfit(data=sbpfit_results, config=config)

    # output results
    csv_writer = csv.writer(sys.stdout)
    csv_writer.writerow(results.keys())
    csv_writer.writerow(results.values())


if __name__ == "__main__":
    main()
