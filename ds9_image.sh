#!/bin/sh
#
# Copyright (c) 2016-2017 Aaron LI
# MIT license
#
# Use DS9 to visualize the ACIS image, with useful/handy
# arguments/options passed.
# Also touch the output image filename for quicker save.
#
# NOTE:
# Tool `manifest.py` is part of the `chandra-acis-analysis` repository.
#
# Created: 2016-04-16
#
# Change logs:
# 2017-02-13:
#   * Use `manifest.py` from `chandra-acis-analysis`
#   * Some cleanups
#

# Default parameters
DS9_GEOMETRY=${DS9_GEOMETRY:-1288x1026-0+0}
DS9_REG_FORMAT=${DS9_REG_FORMAT:-ds9}
DS9_SMOOTH_RADIUS=${DS9_SMOOTH_RADIUS:-4}
DS9_SCALE=${DS9_SCALE:-asinh}
DS9_CMAP=${DS9_CMAP:-sls}


case "$1" in
    -[hH]*)
        echo "Usage:"
        echo "    `basename $0` <img_dir1> ..."
        exit 1
        ;;
esac


INIT_DIR=$(pwd -P)
while [ ! -z "$1" ]; do
    imgdir="$1"
    shift
    cd ${imgdir}
    echo "====== ${PWD} ======"
    FILE=$(manifest.py getpath img_fill)
    if [ -z "${FILE}" ]; then
        echo "*** WARNING: no image file found ***"
        continue
    fi
    echo "FITS file: ${FILE}"
    #
    PNG_FILE="$(echo "${FILE%.fits}" | sed -e 's/_c7//' -e 's/_c0-3//').png"
    [ ! -f "${PNG_FILE}" ] && touch ${PNG_FILE}
    echo "PNG file: ${PNG_FILE}"
    # FoV region
    FOV=$(manifest.py getpath fov)
    DS9_REG_FOV="-regions color green -regions ${REG_FOV}"
    #
    ds9 ${FILE} \
        -bin factor 1 \
        -smooth radius ${DS9_SMOOTH_RADIUS} -smooth yes \
        -scale ${DS9_SCALE} \
        -cmap ${DS9_CMAP} \
        -geometry ${DS9_GEOMETRY} \
        -regions format ${DS9_REG_FORMAT} \
        -regions color white -regions ${REG_FILE} \
        ${DS9_REG_FOV}
    #
    cd ${INIT_DIR}
done
