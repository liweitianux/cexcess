#!/usr/bin/env python3
#
# Copyright (c) 2017 Aaron LI
# MIT license

"""
Get the SDSS image around a specified position using the
"Image Cutout" service.

* Image Cutout:
  http://skyserver.sdss.org/dr13/en/help/docs/api.aspx#imgcutout
* Image List Tool help:
  http://skyserver.sdss.org/dr13/en/tools/started/list.aspx
"""

import os
import sys
import argparse

import requests


def get_image(outfile, ra, dec, scale=0.4, width=1024, height=1024,
              opt="", clobber=False):
    """
    Get and save the requested SDSS image.

    Credits:
    [1] Requests - Raw Response Content
        http://docs.python-requests.org/en/master/user/quickstart/
    [2] How to download image using requests?
        https://stackoverflow.com/a/13137873/4856091
    """
    if os.path.exists(outfile) and clobber is False:
        raise OSError("output file already exists: %s" % outfile)

    params = {
        "ra": ra,
        "dec": dec,
        "scale": scale,
        "width": width,
        "height": height,
        "opt": opt
    }
    url = "http://skyserver.sdss.org/dr13/SkyServerWS/ImgCutout/getjpeg"
    print("Get image from SDSS ...", file=sys.stderr)
    r = requests.get(url, params=params, stream=True)
    print("URL:", r.url, file=sys.stderr)
    r.raise_for_status()
    with open(outfile, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024):
            f.write(chunk)
    print("Saved image to file:", outfile, file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Get SDSS image")
    parser.add_argument("-o", "--outfile", dest="outfile",
                        help="output image filename [jpeg]")
    parser.add_argument("-C", "--clobber", dest="clobber", action="store_true",
                        help="overwrite existing file")
    parser.add_argument("-r", "--ra", dest="ra",
                        type=float, required=True,
                        help="right ascention [degree]")
    parser.add_argument("-d", "--dec", dest="dec",
                        type=float, required=True,
                        help="declination [degree]")
    parser.add_argument("-s", "--scale", dest="scale",
                        type=float, default=0.4,
                        help="scale of image [arcsec/pixel] (default: 0.4)")
    parser.add_argument("-W", "--width", dest="width",
                        type=int, default=1024,
                        help="image width [pixel] (default: 1024)")
    parser.add_argument("-H", "--height", dest="height",
                        type=int, default=1024,
                        help="image height [pixel] (default: 1024)")
    parser.add_argument("-G", "--opt-G", dest="optG", action="store_true",
                        help="overlay coordinate grid and scale")
    parser.add_argument("-L", "--opt-L", dest="optL", action="store_true",
                        help="overlay label of RA/Dec, scale and zoom factor")
    parser.add_argument("-P", "--opt-P", dest="optP", action="store_true",
                        help="mark photometrically identified objects " +
                        "with blue circles")
    parser.add_argument("-S", "--opt-S", dest="optS", action="store_true",
                        help="mark objects with measured spectra " +
                        "using red squares")
    parser.add_argument("-O", "--opt-O", dest="optO", action="store_true",
                        help="mark outlines of photometrically " +
                        "identified objects")
    parser.add_argument("-B", "--opt-B", dest="optB", action="store_true",
                        help="put a box around the outline for each object")
    parser.add_argument("-F", "--opt-F", dest="optF", action="store_true",
                        help="mark boundaries of SDSS fields")
    parser.add_argument("-M", "--opt-M", dest="optM", action="store_true",
                        help="mark regions masked out of the science data")
    parser.add_argument("-Q", "--opt-Q", dest="optQ", action="store_true",
                        help="show boundaries of SDSS spectroscopic plates")
    parser.add_argument("-I", "--opt-I", dest="optI", action="store_true",
                        help="invert the colors (i.e., white background)")
    args = parser.parse_args()

    options = ""
    for optv, optc in [(args.optG, "G"),
                       (args.optL, "L"),
                       (args.optP, "P"),
                       (args.optS, "S"),
                       (args.optO, "O"),
                       (args.optB, "B"),
                       (args.optF, "F"),
                       (args.optM, "M"),
                       (args.optQ, "Q"),
                       (args.optI, "I")]:
        if optv:
            options += optc

    if args.outfile:
        outfile = args.outfile
    else:
        outfile = "sdss_ra{ra:.2f}_dec{dec:.2f}_s{scale:.1f}_{opt}.jpg".format(
            ra=args.ra, dec=args.dec, scale=args.scale, opt=options)

    print("RA:", args.ra, file=sys.stderr)
    print("Dec:", args.dec, file=sys.stderr)
    print("Scale:", args.scale, file=sys.stderr)
    print("Size: %dx%d" % (args.width, args.height), file=sys.stderr)
    print("Options:", options, file=sys.stderr)
    print("Outfile:", outfile, file=sys.stderr)

    get_image(outfile=outfile, ra=args.ra, dec=args.dec, scale=args.scale,
              width=args.width, height=args.height, opt=options,
              clobber=args.clobber)


if __name__ == "__main__":
    main()
