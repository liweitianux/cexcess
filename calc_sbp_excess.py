#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Aaron LI
# Created: 2016-04-26
# Updated: 2016-04-26
#
# TODO:
#   * to estimate errors for the excess value and ratio (Monte Carlo:
#     repeatedly fit the randomized SBP data and calculate excess)
#

"""
Calculate the central brightness excess value and ratio with respect to the
fitted SBP model (i.e., single-beta model).

NOTE:
* excess value: brightness_observed - brightness_model_predicted
* excess ratio: excess_value / brightness_observed
"""

__version__ = "0.1.0"
__date__    = "2016-04-26"


import sys
import argparse
import json
from collections import OrderedDict

import numpy as np
from configobj import ConfigObj

from sbp_fit import make_model


def calc_excess(data, model):
    """
    Calculate the central brightness excess value and ratio with respect
    to the fitted SBP model (single-beta).

    Arguments:
      * data: raw 4-column SBP data
      * model: fitted SBP model
    """
    radii            = data[:, 0]
    radii_err        = data[:, 1]
    brightness       = data[:, 2]
    brightness_model = model.f(radii)
    rin   = radii - radii_err
    rout  = radii + radii_err
    areas = np.pi * (rout**2 - rin**2)
    bsum_obs   = np.sum(brightness * areas)
    bsum_model = np.sum(brightness_model * areas)
    excess_value = bsum_obs - bsum_model
    excess_ratio = excess_value / bsum_obs
    excess = {
            "brightness_obs":   bsum_obs,
            "brightness_model": bsum_model,
            "excess_value":     excess_value,
            "excess_ratio":     excess_ratio,
    }
    return excess


def main():
    # default arguments
    default_outfile = "excess.json"

    parser = argparse.ArgumentParser(
            description="Calculate the central brightness excess value",
            epilog="Version: %s (%s)" % (__version__, __date__))
    parser.add_argument("-V", "--version", action="version",
            version="%(prog)s " + "%s (%s)" % (__version__, __date__))
    parser.add_argument("config", help="Config file for SBP fitting")
    parser.add_argument("outfile", nargs="?", default=default_outfile,
            help="Output json file to save the results " + \
                 "(default: %s)" % default_outfile)
    args = parser.parse_args()

    config = ConfigObj(args.config)
    modelname = config["model"]
    config_model = config[modelname]

    model = make_model(config, modelname=modelname)
    print("SBP model: %s" % model.name, file=sys.stderr)

    sbpfit_outfile = config.get("outfile")
    sbpfit_outfile = config_model.get("outfile", sbpfit_outfile)
    sbpfit_results = json.load(open(sbpfit_outfile),
                               object_pairs_hook=OrderedDict)

    # Load fitted parameters for model
    for pname, pvalue in sbpfit_results["params"].items():
        model.set_param(name=pname, value=pvalue[0])

    sbpdata = np.loadtxt(config["sbpfile"])
    excess = calc_excess(data=sbpdata, model=model)

    excess_data = OrderedDict([
            ("name",             config["name"]),
            ("obsid",            int(config["obsid"])),
            ("model",            modelname),
            ("brightness_obs",   excess["brightness_obs"]),
            ("brightness_model", excess["brightness_model"]),
            ("excess_value",     excess["excess_value"]),
            ("excess_ratio",     excess["excess_ratio"]),
    ])
    excess_data_json = json.dumps(excess_data, indent=2)
    print(excess_data_json)
    open(args.outfile, "w").write(excess_data_json+"\n")


if __name__ == "__main__":
    main()

#  vim: set ts=4 sw=4 tw=0 fenc=utf-8 ft=python: #
