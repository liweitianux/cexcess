#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Aaron LI
# Created: 2016-04-26
# Updated: 2016-05-17
#
# Change log:
# 2016-05-17:
#   * Add argument "--subtract-bkg" and consider background subtraction
#   * Add argument "--r500-cut" and `rcut` support
# 2016-05-06:
#   * Add function `estimate_excess_error()` to estimate uncertainty
#   * update according to `sbp_fit` renamed to `fit_sbp`
#   * PEP8 fixes
#

"""
Calculate the central brightness excess value and ratio with respect to the
fitted SBP model (i.e., single-beta model).

NOTE:
* excess value: brightness_observed - brightness_model_predicted
* excess ratio: excess_value / brightness_observed
"""


import sys
import argparse
import json
from collections import OrderedDict

import numpy as np
from configobj import ConfigObj

from fit_sbp import make_model, make_sbpfit

__version__ = "0.3.1"
__date__ = "2016-05-17"


def calc_excess(data, fitted_model, rcut=None,
                subtract_bkg=False, verbose=False):
    """
    Calculate the central brightness excess value and ratio with respect
    to the fitted SBP model (single-beta).

    TODO/XXX:
      * whether to interpolate the SBP?

    Arguments:
      * data: raw 4-column SBP data
      * fitted_model: fitted SBP model
      * rcut: cut radius for total/integrated brightness calculation;
              if rcut larger than the maximum available radius, then
              use the maximum radius instead.
      * subtract_bkg: whether subtract the fitted background?
    """
    radii = data[:, 0]
    radii_err = data[:, 1]
    brightness = data[:, 2]
    brightness_model = fitted_model.f(radii)
    rin = radii - radii_err
    rout = radii + radii_err
    if rcut is not None and rcut < rout[-1]:
        ncut = np.sum(rin <= rcut)
        rin = rin[:ncut]
        rout = rout[:ncut]
        rout[-1] = rcut
        brightness = brightness[:ncut]
        brightness_model = brightness_model[:ncut]
        if verbose:
            print("DEBUG: rcut:", rcut, file=sys.stderr)
            print("DEBUG: ncut:", ncut, file=sys.stderr)
            print("DEBUG: rin:", rin, file=sys.stderr)
            print("DEBUG: rout:", rout, file=sys.stderr)
            print("DEBUG: brightness:", brightness, file=sys.stderr)
    if subtract_bkg:
        bkg = fitted_model.get_param("bkg").value
        if verbose:
            print("Subtract fitted background: %g" % bkg)
        brightness -= bkg
        brightness_model -= bkg
    areas = np.pi * (rout**2 - rin**2)
    bsum_obs = np.sum(brightness * areas)
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


def estimate_excess_error(data, sbpfit, mctimes, rcut=None,
                          subtract_bkg=False, verbose=False):
    """
    Estimate the uncertainty of central excess value by Monte Carlo method.

    XXX/TODO:
      * whether also consider the uncertainty of R500?

    Arguments:
      * data: 4-column SBP data (radius, r_err, brightness, brightness_err)
      * sbpfit: `SbpFit` object used to perform SBP fitting
      * mctimes: number of Monte Carlo iterations
    """
    brightness = data[:, 2]
    brightness_err = data[:, 3]
    params = sbpfit.dump_params()
    ev_results = []
    er_results = []
    if verbose:
        print("Estimating excess uncertainty by Monte Carlo " +
              "(%d times) ..." % mctimes, end="", flush=True)
    for i in range(mctimes):
        if verbose and i % 100 == 0:
            print("%d..." % i, end="", flush=True)
        # randomize SBP data
        brightness_rand = np.random.normal(brightness, scale=brightness_err)
        sbpdata_rand = data.copy()
        sbpdata_rand[:, 2] = brightness_rand
        # load randomized data and perform SBP fit
        sbpfit.reset(keep_source=True)
        sbpfit.load_data(sbpdata_rand, keep_mask=True)
        sbpfit.load_params(params)
        sbpfit.fit()
        model_rand = sbpfit.get_model()
        excess = calc_excess(data=sbpdata_rand, fitted_model=model_rand,
                             rcut=rcut, subtract_bkg=subtract_bkg,
                             verbose=False)
        ev_results.append(excess["excess_value"])
        er_results.append(excess["excess_ratio"])
    if verbose:
        print("DONE!", flush=True)
    ev_mean = np.mean(ev_results)
    ev_std = np.std(ev_results)
    ev_q25, ev_median, ev_q75 = np.percentile(ev_results, q=(25, 50, 75))
    er_mean = np.mean(er_results)
    er_std = np.std(er_results)
    er_q25, er_median, er_q75 = np.percentile(er_results, q=(25, 50, 75))
    results = {
        "excess_value_mean":   ev_mean,
        "excess_value_median": ev_median,
        "excess_value_q25":    ev_q25,
        "excess_value_q75":    ev_q75,
        "excess_value_std":    ev_std,
        "excess_ratio_mean":   er_mean,
        "excess_ratio_median": er_median,
        "excess_ratio_q25":    er_q25,
        "excess_ratio_q75":    er_q75,
        "excess_ratio_std":    er_std,
    }
    return results


def main():
    # default arguments
    default_outfile = "excess.json"
    default_mctimes = 1000

    parser = argparse.ArgumentParser(
            description="Calculate the central brightness excess value",
            epilog="Version: %s (%s)" % (__version__, __date__))
    parser.add_argument("-V", "--version", action="version",
                        version="%(prog)s " +
                                "%s (%s)" % (__version__, __date__))
    parser.add_argument("-v", "--verbose", dest="verbose",
                        required=False, action="store_true",
                        help="show verbose information")
    parser.add_argument("-m", "--mctimes", dest="mctimes", required=False,
                        type=int, default=default_mctimes,
                        help="number of MC times to estimate excess error " +
                             "(default: %d)" % default_mctimes)
    parser.add_argument("-R", "--r500-cut", dest="r500_cut",
                        type=float, default=0.5,
                        help="fraction of R500 to be taken as the cut " +
                             "radius for total brightness calculation " +
                             "(default: 0.5)")
    parser.add_argument("-B", "--subtract-bkg", dest="subtract_bkg",
                        action="store_true",
                        help="subtract the fitted background and calculate " +
                             "the net brightness")
    parser.add_argument("config", help="Config file for SBP fitting")
    parser.add_argument("outfile", nargs="?", default=default_outfile,
                        help="Output json file to save the results " +
                             "(default: %s)" % default_outfile)
    args = parser.parse_args()

    config = ConfigObj(args.config)

    r500_pix = float(config["r500_pix"])
    rcut = args.r500_cut * r500_pix

    modelname = config["model"]
    config_model = config[modelname]

    model = make_model(config, modelname=modelname)
    print("SBP model: %s" % model.long_name, file=sys.stderr)

    sbpfit_outfile = config.get("outfile")
    sbpfit_outfile = config_model.get("outfile", sbpfit_outfile)
    sbpfit_results = json.load(open(sbpfit_outfile),
                               object_pairs_hook=OrderedDict)

    # Load fitted parameters for model
    for pname, pvalue in sbpfit_results["params"].items():
        model.set_param(name=pname, value=pvalue[0])

    sbpfit = make_sbpfit(config, model=model)

    sbpdata = np.loadtxt(config["sbpfile"])
    excess = calc_excess(data=sbpdata, fitted_model=model, rcut=rcut,
                         subtract_bkg=args.subtract_bkg,
                         verbose=args.verbose)
    excess_err = estimate_excess_error(data=sbpdata, sbpfit=sbpfit,
                                       mctimes=args.mctimes, rcut=rcut,
                                       subtract_bkg=args.subtract_bkg,
                                       verbose=args.verbose)

    excess_data = OrderedDict([
        ("name",                config["name"]),
        ("obsid",               int(config["obsid"])),
        ("model",               modelname),
        ("excess_rcut",         rcut),
        ("brightness_obs",      excess["brightness_obs"]),
        ("brightness_model",    excess["brightness_model"]),
        ("excess_value",        excess["excess_value"]),
        ("excess_value_mean",   excess_err["excess_value_mean"]),
        ("excess_value_median", excess_err["excess_value_median"]),
        ("excess_value_q25",    excess_err["excess_value_q25"]),
        ("excess_value_q75",    excess_err["excess_value_q75"]),
        ("excess_value_std",    excess_err["excess_value_std"]),
        ("excess_ratio",        excess["excess_ratio"]),
        ("excess_ratio_mean",   excess_err["excess_ratio_mean"]),
        ("excess_ratio_median", excess_err["excess_ratio_median"]),
        ("excess_ratio_q25",    excess_err["excess_ratio_q25"]),
        ("excess_ratio_q75",    excess_err["excess_ratio_q75"]),
        ("excess_ratio_std",    excess_err["excess_ratio_std"]),
        ("mc_times",            args.mctimes),
    ])
    excess_data_json = json.dumps(excess_data, indent=2)
    print(excess_data_json)
    open(args.outfile, "w").write(excess_data_json+"\n")


if __name__ == "__main__":
    main()

#  vim: set ts=4 sw=4 tw=0 fenc=utf-8 ft=python: #
