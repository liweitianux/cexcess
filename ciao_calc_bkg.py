#!/usr/bin/env python3
#
# Calculate the *background surface brightness (SB)* level from the
# *corrected background spectrum*.
# The calculated background SB value is used to provide constraint for
# surface brightness profile (SBP) fitting, and is also adopted to
# subtract the background contribution before carrying out the SBP
# deprojection.
#
# Aaron LI
# Created: 2016-06-10
# Updated: 2016-06-10
#

import argparse
import subprocess
import tempfile
import json
from collections import OrderedDict

import numpy as np
from astropy.io import fits


def energy2channel(energy):
    """
    Convert energy (eV) to Chandra ACIS channel number.

    Reference: http://cxc.harvard.edu/ciao/dictionary/pi.html
    """
    return int((energy / 14.6) + 1)


def parse_erange(erange):
    """
    Parse the given erange string, and return the lower and upper energies.
    """
    e_low, e_high = map(float, erange.split(":"))
    return (e_low, e_high)


def calc_exp(expmap, regfile):
    """
    Calculate the area of the background spectrum extraction region,
    and the mean exposure (un-normalized w.r.t exposure time) of that region.
    """
    tf = tempfile.NamedTemporaryFile()
    cmd_args = [
        "dmextract",
        "infile=%s[bin sky=region(%s)]" % (expmap, regfile),
        "outfile=%s" % tf.name,
        "opt=generic", "clobber=yes"
    ]
    subprocess.run(["punlearn", "dmextract"])
    subprocess.run(args=cmd_args)
    with fits.open(tf.name) as hist_fits:
        total_exp = hist_fits["HISTOGRAM"].data["COUNTS"][0]
        total_exp_err = hist_fits["HISTOGRAM"].data["ERR_COUNTS"][0]
        area = hist_fits["HISTOGRAM"].data["AREA"][0]
        mean_exp = hist_fits["HISTOGRAM"].data["SUR_BRI"][0]
        mean_exp_err = hist_fits["HISTOGRAM"].data["SUR_BRI_ERR"][0]
    tf.close()
    return {
        "total_exp":     total_exp,
        "total_exp_err": total_exp_err,
        "area":          area,
        "mean_exp":      mean_exp,
        "mean_exp_err":  mean_exp_err,
    }


def calc_spec_counts(spec, erange):
    """
    Calculate the spectrum total counts within the specified energy range.
    Also extract the EXPOSURE and BACKSCAL information.
    """
    specfits = fits.open(spec)
    channel = specfits["SPECTRUM"].data["CHANNEL"]
    counts = specfits["SPECTRUM"].data["COUNTS"]
    channel_low, channel_high = map(energy2channel, erange)
    chan_idx = np.where(np.logical_and(channel >= channel_low,
                                       channel <= channel_high))
    total_counts = np.sum(counts[chan_idx])
    return {
        "energy_low": erange[0],
        "energy_high": erange[1],
        "channel_low": channel_low,
        "channel_high": channel_high,
        "counts": int(total_counts),
        "exposure": specfits["SPECTRUM"].header["EXPOSURE"],
        "backscal": specfits["SPECTRUM"].header["BACKSCAL"],
    }


def main():
    parser = argparse.ArgumentParser(
            description="Calculate the background surface brightness")
    parser.add_argument("-b", "--orig-bkg", dest="orig_bkg", required=True,
                        help="original/uncorrected local background " +
                        "spectrum (e.g., lbkg.pi), from which to get the " +
                        "original/source EXPOSURE and BACKSCAL, etc.")
    parser.add_argument("-B", "--corr-bkg", dest="corr_bkg", required=True,
                        help="corrected background spectrum for the " +
                        "Galactic and cosmic background radiations " +
                        "(e.g., bkgcorr_blanksky_lbkg.pi)")
    parser.add_argument("-r", "--bkg-region", dest="bkg_region", required=True,
                        help="region used to extract the background " +
                        "spectrum (e.g., lbkg.reg")
    parser.add_argument("-e", "--expmap", dest="expmap", required=True,
                        help="exposure map w.r.t the spectrum")
    parser.add_argument("-E", "--erange", dest="erange", required=True,
                        help="energy range used for the exposure map " +
                        "(e.g., 700:7000)")
    parser.add_argument("-o", "--outfile", dest="outfile",
                        required=False, default="sb_bkg.json",
                        help="json file to save the background SB results")
    args = parser.parse_args()

    e_low, e_high = parse_erange(args.erange)
    orig_bkg_results = calc_spec_counts(args.orig_bkg, erange=(e_low, e_high))
    corr_bkg_results = calc_spec_counts(args.corr_bkg, erange=(e_low, e_high))
    exp_results = calc_exp(args.expmap, regfile=args.bkg_region)

    corr_counts = corr_bkg_results["counts"]
    corr_exposure = corr_bkg_results["exposure"]
    corr_backscal = corr_bkg_results["backscal"]
    orig_exposure = orig_bkg_results["exposure"]
    orig_backscal = orig_bkg_results["backscal"]
    area = exp_results["area"]
    mean_exp = exp_results["mean_exp"]
    bkg_sb = corr_counts / corr_exposure / (mean_exp / orig_exposure) \
             / (area * corr_backscal / orig_backscal)

    results = OrderedDict([
        ("energy_low",    e_low),
        ("energy_high",   e_high),
        ("channel_low",   corr_bkg_results["channel_low"]),
        ("channel_high",  corr_bkg_results["channel_high"]),
        ("counts",        corr_counts),
        ("exposure",      orig_exposure),
        ("exposure_bkg",  corr_exposure),
        ("backscal",      corr_backscal),
        ("total_exp",     exp_results["total_exp"]),
        ("total_exp_err", exp_results["total_exp_err"]),
        ("area",          exp_results["area"]),
        ("mean_exp",      exp_results["mean_exp"]),
        ("mean_exp_err",  exp_results["mean_exp_err"]),
        ("bkg_sb",        bkg_sb),
    ])
    results_json = json.dumps(results, indent=2)
    print(results_json)
    open(args.outfile, "w").write(results_json+"\n")


if __name__ == "__main__":
    main()
