#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Aaron LI
# Created: 2016-04-29
# Updated: 2016-05-02
#
# TODO:
#   * to calculate the PEI error
#

"""
Calculate the power excess index (PEI), which is defined the area ratio of
the lower-left part with respect to the total rectangle, which is further
defined by the radial power spectrum and the scale of 0.035R500 and 0.35R500,
in the logarithmic space.

Reference:
Zhang, C., et al. 2016, ApJ
"""


import os
import glob
import argparse
import json
from collections import OrderedDict

import numpy as np
import scipy.interpolate
import scipy.integrate

import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.path import Path
import matplotlib.patches as patches

from make_r500_regions import get_r500

__version__ = "0.3.1"
__date__ = "2016-05-02"

plt.style.use("ggplot")


def calc_pei(data, r500, interp_np=101):
    """
    Calculate the power excess index (PEI), which is defined the area ratio
    of the lower-left part with respect to the total rectangle.

    Arguments:
      * data: 3-column power spectrum data (frequency, power, power_err)
      * r500: R500 value in unit of the inverse of the above "frequency"
      * interp_np: number of data points interpolated to calculate PEI
    """
    freqs = data[:, 0]
    psd1d = data[:, 1]
    # frequency values corresponding to 0.35R500 and 0.035R500
    pei_fmin = 1.0 / (0.350 * r500)
    pei_fmax = 1.0 / (0.035 * r500)
    # switch to the logarithmic scale
    # XXX: how to deal with the errors
    mask = (freqs > 0.0)
    x = np.log10(freqs[mask])
    y = np.log10(psd1d[mask])
    pei_xmin = np.log10(pei_fmin)
    pei_xmax = np.log10(pei_fmax)
    # interpolate the power spectrum
    f_interp = scipy.interpolate.interp1d(x, y, kind="linear",
                                          assume_sorted=True)
    x_interp = np.linspace(pei_xmin, pei_xmax, num=interp_np)
    y_interp = f_interp(x_interp)
    pei_ymin = np.min(y_interp)
    pei_ymax = np.max(y_interp)
    # calculate the PEI
    area_total = (pei_xmax - pei_xmin) * (pei_ymax - pei_ymin)
    area_below = scipy.integrate.trapz((y_interp-pei_ymin), x_interp)
    pei_value = area_below / area_total
    results = {
            "area_total": area_total,
            "area_below": area_below,
            "pei_value":  pei_value,
            "pei_err":    None,
    }
    data_interp_log10 = np.column_stack((x_interp, y_interp))
    return (results, data_interp_log10)


def plot_pei(data, data_interp_log10, info={}, ax=None, fig=None):
    """
    Make a plot to visualize the PEI rectangular.
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    # prepare data
    freqs = data[:, 0]
    psd1d = data[:, 1]
    psd1d_err = data[:, 2]
    x_interp = 10 ** data_interp_log10[:, 0]
    y_interp = 10 ** data_interp_log10[:, 1]
    pei_xmin = np.min(x_interp)
    pei_xmax = np.max(x_interp)
    pei_ymin = np.min(y_interp)
    pei_ymax = np.max(y_interp)
    #
    mask = (freqs > 0.0)
    xmin = np.min(freqs[mask]) / 1.2
    xmax = np.max(freqs[mask])
    ymin = np.min(psd1d[mask]) / 3.0
    ymax = np.max(psd1d[mask] + psd1d_err[mask]) * 1.2
    #
    ax.plot(freqs, psd1d, color="black", linestyle="none",
            marker="o", markersize=5, alpha=0.7)
    ax.errorbar(freqs, psd1d, yerr=psd1d_err, fmt="none",
                ecolor="blue", alpha=0.7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title("Power Spectral Density & PEI (%s; %d)" %
                 (info.get("name"), info.get("obsid")))
    ax.set_xlabel(r"k (pixel$^{-1}$)")
    ax.set_ylabel("Power")
    ax.text(x=xmax/1.1, y=ymax/1.2,
            s="PEI = %.2f / %.2f = %.2f" % (info.get("area_below"),
                                            info.get("area_total"),
                                            info.get("pei")),
            horizontalalignment="right", verticalalignment="top")
    # plot the interpolated data points and the PEI rectangle
    # credit: http://matplotlib.org/users/path_tutorial.html
    ax.plot(x_interp, y_interp, linestyle="--", marker="D", markersize=2,
            color="green", alpha=0.9)
    vertices = [
        (pei_xmin, pei_ymin),  # left, bottom
        (pei_xmin, pei_ymax),  # left, top
        (pei_xmax, pei_ymax),  # right, top
        (pei_xmax, pei_ymin),  # right, bottom
        (pei_xmin, pei_ymin),  # ignored
    ]
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    path = Path(vertices, codes)
    patch = patches.PathPatch(path, fill=False, color="green", linewidth=2,
                              alpha=0.9)
    ax.add_patch(patch)
    fig.tight_layout()
    return (fig, ax)


def main():
    # default arguments
    default_infile = "psd.txt"
    default_outfile = "pei.json"
    default_infojson = "../*_INFO.json"

    parser = argparse.ArgumentParser(
            description="Calculate the power excess index (PEI)",
            epilog="Version: %s (%s)" % (__version__, __date__))
    parser.add_argument("-V", "--version", action="version",
                        version="%(prog)s " + "%s (%s)" % (__version__,
                                                           __date__))
    parser.add_argument("-j", "--json", dest="json", required=False,
                        help="the *_INFO.json file " +
                             "(default: find %s)" % default_infojson)
    parser.add_argument("-i", "--infile", dest="infile",
                        required=False, default=default_infile,
                        help="input data of the radial power spectrum " +
                             "(default: %s)" % default_infile)
    parser.add_argument("-o", "--outfile", dest="outfile",
                        required=False, default=default_outfile,
                        help="output json file to save the results " +
                             "(default: %s)" % default_outfile)
    parser.add_argument("-p", "--png", dest="png", default=None,
                        help="make PEI plot and save " +
                             "(default: same basename as outfile)")
    args = parser.parse_args()

    if args.png is None:
        args.png = os.path.splitext(args.outfile)[0] + ".png"

    info_json = glob.glob(default_infojson)[0]
    if args.json:
        info_json = args.json

    json_str = open(info_json).read().rstrip().rstrip(",")
    info = json.loads(json_str)
    name = info["Source Name"]
    obsid = int(info["Obs. ID"])
    r500 = get_r500(info)
    r500_kpc = r500["r500_kpc"]
    r500_pix = r500["r500_pix"]
    kpc_per_pix = r500["kpc_per_pix"]

    psd_data = np.loadtxt(args.infile)
    pei, data_interp_log10 = calc_pei(psd_data, r500=r500_pix)

    pei_data = OrderedDict([
            ("name",            name),
            ("obsid",           obsid),
            ("r500_kpc",        r500_kpc),
            ("r500_pix",        r500_pix),
            ("kpc_per_pix",     kpc_per_pix),
            ("area_total",      pei["area_total"]),
            ("area_below",      pei["area_below"]),
            ("pei",             pei["pei_value"]),
            ("pei_err",         pei["pei_err"]),
    ])
    pei_data_json = json.dumps(pei_data, indent=2)
    print(pei_data_json)
    open(args.outfile, "w").write(pei_data_json+"\n")

    # Make and save a plot
    fig = Figure(figsize=(10, 8))
    FigureCanvas(fig)
    ax = fig.add_subplot(111)
    plot_pei(psd_data, data_interp_log10, info=pei_data, ax=ax, fig=fig)
    fig.savefig(args.png, format="png", dpi=150)


if __name__ == "__main__":
    main()

#  vim: set ts=4 sw=4 tw=0 fenc=utf-8 ft=python: #
