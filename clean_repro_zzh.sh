#!/bin/sh
#
# Clean the 'repro' and its subdirectories for the sources
# drawn from zzh among the excess_sample.
#
#
# Aaron LI
# 2016-04-12
#

case "$1" in
    ""|-[hH]*)
        echo "Usage: `basename $0` <repro_dir1> ..."
        exit 1
esac

INIT_DIR=`pwd -P`
while [ ! -z "$1" ]; do
    repro_dir="$1"
    shift
    cd ${INIT_DIR}
    cd ${repro_dir}
    echo "*** Cleaning '${repro_dir}' ..."
    cd evt
    rm -fv *~ .?*~ _* *.log
    rm -fv evt2*_orig.fits evt2*_rmsrcs.fits
    rm -fv img_*.fits sbprofile.reg *fov*.fits
    mv -fv peak*.reg ../img/
    cd ../bkg
    rm -fv *~ .?*~ _* *.log
    cd ../img
    rm -fv *~ .?*~ _* *.log
    rm -fv rspec*bak *fov*.fits img_*.fits
done

