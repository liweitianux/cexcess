#!/usr/bin/env python3
#
# Extract surface brightness profile, with optional background
# subtraction and region exclusion.
#
# Aaron LI
# Created: 2016-05-17
#
# Change log:
#

import os
import re
import argparse
import subprocess
import tempfile

import numpy as np
from astropy.io import fits


def check_acis_type(infile):
    """
    Check the ACIS type of the infile: ACIS-I or ACIS-S
    """
    subprocess.run(["punlearn", "dmkeypar"])
    ret = subprocess.run(["dmkeypar", infile, "DETNAM", "echo=yes"],
                         check=True, stdout=subprocess.PIPE)
    detnam = ret.stdout.decode("utf-8")
    if re.match(r"^ACIS-0123[4-9]*", detnam):
        results = {"type": "ACIS-I", "ccd": "0:3"}
    elif re.match(r"^ACIS-[0-6]*7", detnam):
        results = {"type": "ACIS-S", "ccd": "7"}
    else:
        raise ValueError("unknown DETNAM: %s" % detnam)
    return results


def make_sbp_region(regfile, exclude_regfile=None, fov=None, ccd=None):
    """
    Make the regions for SBP extraction, considering the regions to be
    excluded and FoV constraint.

    Return a list containing all the SBP regions.
    """
    regions = map(str.strip, open(regfile).readlines())
    regions = list(filter(lambda x: re.match(r"^(circle|annulus|pie).*",
                                             x, re.I),
                          regions))
    if exclude_regfile is not None:
        ex_regions = map(str.strip, open(exclude_regfile).readlines())
        ex_regions = list(filter(lambda x: not re.match(r"^\s*(|#.*)\s*$", x),
                                 ex_regions))
        ex_regions = "*!".join([""] + ex_regions)
        regions = [reg+ex_regions for reg in regions]
    if fov is not None:
        with tempfile.NamedTemporaryFile() as fp:
            subprocess.run(["punlearn", "dmmakereg"])
            subprocess.run(args=[
                "dmmakereg",
                "region=region(%s[ccd_id=%s])" % (fov, ccd),
                "outfile=%s" % fp.name,
                "kernel=ascii",
                "clobber=yes"
            ])
            fov_regions = map(str.strip, open(fp.name).readlines())
        fov_regions = list(filter(lambda x: re.match(r"^physical;polygon.*$",
                                                     x, re.I),
                                  fov_regions))
        fov_regions = [re.sub(r"^physical;\s*", "",
                              re.sub(r"\s*#\s*$", "", reg))
                       for reg in fov_regions]
        regions = [
            "+".join([
                reg + "*" + fov for fov in fov_regions
            ])
            for reg in regions
        ]
    return regions


def extract_sbp(infile, expmap, regfile, outprefix, bkg=None, erange=None):
    """
    Extract the surface brightness profile

    If `bkg` is provided, then background subtraction is considered.
    """
    if erange is not None:
        erange = "[energy=%s]" % erange
    else:
        erange = ""
    sbpfile = outprefix + ".fits"
    cmd_args = [
        "dmextract",
        "infile=%s%s[bin sky=@%s]" % (infile, erange, regfile),
        "outfile=%s" % sbpfile,
        "exp=%s" % expmap,
        "opt=generic",
        "clobber=yes"
    ]
    if bkg is not None:
        # consider background subtraction
        subprocess.run(["punlearn", "dmkeypar"])
        ret = subprocess.run(args=["dmkeypar", infile, "EXPOSURE", "echo=yes"],
                             check=True, stdout=subprocess.PIPE)
        exposure_evt = float(ret.stdout.decode("utf-8"))
        ret = subprocess.run(args=["dmkeypar", bkg, "EXPOSURE", "echo=yes"],
                             check=True, stdout=subprocess.PIPE)
        exposure_bkg = float(ret.stdout.decode("utf-8"))
        bkg_norm = exposure_evt / exposure_bkg
        cmd_args += [
            "bkg=%s%s[bin sky=@%s]" % (bkg, erange, regfile),
            "bkgnorm=%s" % bkg_norm,
            "bkgexp=)exp"
        ]
    subprocess.run(["punlearn", "dmextract"])
    subprocess.run(args=cmd_args)
    # add `RMID` and `R_ERR` columns
    sbpfile_rmid = outprefix + "_rmid.fits"
    subprocess.run(["punlearn", "dmtcalc"])
    subprocess.run(args=[
        "dmtcalc",
        "infile=%s" % sbpfile,
        "outfile=%s" % sbpfile_rmid,
        "expression=RMID=(R[0]+R[1])/2,R_ERR=(R[1]-R[0])/2",
        "clobber=yes"
    ])
    # dump SBP data
    with fits.open(sbpfile_rmid) as sbpfits:
        rmid = sbpfits["HISTOGRAM"].data["RMID"]
        r_err = sbpfits["HISTOGRAM"].data["R_ERR"]
        sur_flux = sbpfits["HISTOGRAM"].data["SUR_FLUX"]
        sur_flux_err = sbpfits["HISTOGRAM"].data["SUR_FLUX_ERR"]
    sbpdata = np.column_stack([rmid, r_err, sur_flux, sur_flux_err])
    sbp_txt = outprefix + ".txt"
    np.savetxt(sbp_txt, sbpdata, header="RMID  R_ERR  SUR_FLUX  SUR_FLUX_ERR")
    # create a QDP file
    sbp_qdp = map(str.strip, open(sbp_txt).readlines())
    sbp_qdp = [re.sub(r"#", "!", line) for line in sbp_qdp]
    sbp_qdp = [
        "READ SERR 1 2",
        'LABEL T "Surface Brightness Profile"',
        'LABEL X "Radius (pixel)"',
        'LABEL Y "Surface Flux (photons/cm\\u2\\d/pixel\\u2\\d/s)"',
        "LOG X Y ON"
    ] + sbp_qdp
    open(outprefix + ".qdp", "w").write("\n".join(sbp_qdp) + "\n")


def main():
    parser = argparse.ArgumentParser(
            description="Extract surface brightness profile")
    parser.add_argument("-j", "--json", dest="json", required=False,
                        help="the *_INFO.json file " +
                             "(default: find ../*_INFO.json)")
    parser.add_argument("-r", "--region", dest="region",
                        required=False, default="sbprofile.reg",
                        help="surface brightness profile region file " +
                             "(default: sbprofile.reg)")
    parser.add_argument("-R", "--exclude-region", dest="exclude_region",
                        default=None,
                        help="file containing regions to be excluded")
    parser.add_argument("-i", "--infile", dest="infile", required=True,
                        help="input EVT2/image file")
    parser.add_argument("-e", "--expmap", dest="expmap", required=True,
                        help="exposure map of the input file")
    parser.add_argument("-b", "--bkg", dest="bkg", default=None,
                        help="background event/image of the input file")
    parser.add_argument("-E", "--erange", dest="erange", default=None,
                        help="energy range of interest (event files only)")
    parser.add_argument("-F", "--fov", dest="fov", default=None,
                        help="FoV FITS file (for applying FoV constraint " +
                             "to the SBP regions)")
    parser.add_argument("-o", "--outprefix", dest="outprefix",
                        required=False,
                        help="prefix of output files (default: same " +
                             "basename as the input region file)")
    args = parser.parse_args()

    # set output prefix if not specified
    if not args.outprefix:
        args.outprefix = os.path.splitext(args.region)[0]

    acis_type = check_acis_type(args.infile)
    regions = make_sbp_region(regfile=args.region,
                              exclude_regfile=args.exclude_region,
                              fov=args.fov, ccd=acis_type["ccd"])
    regfile_out = args.outprefix + "_fix.reg"
    open(regfile_out, "w").write("\n".join(regions) + "\n")

    extract_sbp(infile=args.infile, expmap=args.expmap,
                regfile=regfile_out, outprefix=args.outprefix,
                bkg=args.bkg, erange=args.erange)


if __name__ == "__main__":
    main()
