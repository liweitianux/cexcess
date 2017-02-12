#!/bin/sh
#
# Copyright (c) 2016 Aaron LI
# MIT license
#
# Create the new `sbpfit' subdirectory, and prepare the files for fitting
# the surface brightness profile.
#
# Created: 2016-03-28
#
# Changelog:
# 2016-04-26:
#   * Do NOT prepare the sbpfit config, just link needed files
#   * Remove the case that source has "ne" and "sw" two parts (dropped sources)
#


case "$1" in
    -[hH]*)
        echo "Usage:"
        echo "    `basename $0` <reprodir1> ..."
        exit 1
        ;;
esac

SBPFIT_DIR="sbpfit"
INIT_DIR=`pwd -P`

while [ ! -z "$1" ]; do
    reprodir="$1"
    shift
    echo "*** ${reprodir} ***"
    cd ${INIT_DIR}
    cd ${reprodir}
    [ -d "${SBPFIT_DIR}" ] && mv -fv ${SBPFIT_DIR} ${SBPFIT_DIR}_bak
    mkdir ${SBPFIT_DIR} && cd ${SBPFIT_DIR}
    img_dir="../img"
    ln -sv ${img_dir}/sbprofile.* .
    ln -sv ${img_dir}/sbprofile_rmid.fits .
    ln -sv ${img_dir}/evt2_c*_clean.fits .
    ln -sv ${img_dir}/*_img_*_fill.png .
done
