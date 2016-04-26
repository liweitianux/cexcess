#!/bin/sh
#
# Clean the 'img' directories for the excess_sample.
#
#
# Aaron LI
# 2016-04-11
#

case "$1" in
    ""|-[hH]*)
        echo "Usage: `basename $0` <imgdir1> ..."
        exit 1
esac

INIT_DIR=`pwd -P`
while [ ! -z "$1" ]; do
    imgdir="$1"
    shift
    cd ${INIT_DIR}
    cd ${imgdir}
    echo "*** Cleaning '${imgdir}' ..."
    rm -fv _* *~ .?*~
    rm -fv celld*.reg csb_results*.txt rspec*bak
    rm -fv img_*.fits *fov*.fits ds9.jpg
    rm -fv sbp*.log expcorr*.log
done

