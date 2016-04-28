#!/bin/sh
#
# Clean and merge the two 'img' directories for ZZH's sample data.
#   * 'img': ${name}/${obsid}/evt2/img
#   * 'img2': ${name}/${obsid}/repro/img
# The contents of 'img2' are merged to 'img'.
#
# Aaron LI
# Created: 2016-04-26
#

IMG_DIR="img"
IMG2_DIR="img2"

case "$1" in
    -[hH]*)
        echo "Usage:"
        echo "    `basename $0` <repro_dir1> ..."
        exit 1
        ;;
esac

INIT_DIR=`pwd -P`
while [ ! -z "$1" ]; do
    repro_dir="$1"
    shift
    cd ${INIT_DIR}
    cd ${repro_dir}
    echo "*** ${PWD} ***"
    if [ ! -d "${IMG2_DIR}" ]; then
        echo "WARNING: '${IMG2_DIR}' does not exists; skipped!"
        continue
    fi
    # clean ${IMG_DIR} and ${IMG2_DIR}
    ( cd ${IMG_DIR}; \
        rm -fv _* *.log *bak evt2_*.fits img_*.fits \;
        rm -fv *smooth* *raw* \;
        rm -fv test* tmp* sources* pntsrc* )
    ( cd ${IMG2_DIR}; \
        rm -fv rspec* sbprofile* radius_sbp.txt flux_sbp.txt )
    # merge
    mv -fv ${IMG2_DIR}/* ${IMG_DIR}
    rmdir -v ${IMG2_DIR}
done

