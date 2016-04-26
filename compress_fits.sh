#!/bin/sh
#
# Compress the FITS files that are not used directly.
#
# Aaron LI
# 2016-04-16
#

# whether to compress all big files (useful for dropped sources)
FLAG_ALL="NO"

# compress command
COMPRESS="xz -v"
#COMPRESS="ls -lh"  # test


case "$1" in
    -[hH]*)
        echo "Usage:"
        echo "    `basename $0` [ -a ] <source_dir1> ..."
        exit 1
        ;;
    -[aA]*)
        FLAG_ALL="YES"
        shift
        ;;
esac


while [ ! -z "$1" ]; do
    source="$1"
    shift
    echo "====== ${source} ======"
    find ${source}/ -type f \
        \( -name 'acis*_repro_evt2.fits' -o \
           -name 'pcadf*_asol*.fits' -o \
           -name 'blanksky_*.fits' -o \
           -name '*_tdet.fits' -o \
           -name 'imgcorr_*.fits' \) \
        -exec ${COMPRESS} '{}' \;
    if [ "${FLAG_ALL}" = "YES" ]; then
        echo "*** ALL ***"
        find ${source}/ -type f \
            \( -name 'evt2_*.fits' -o \
               -name 'expmap_*.fits' -o \
               -name 'img_*.fits' \) \
            -exec ${COMPRESS} '{}' \;
    fi
done

