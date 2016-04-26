#!/bin/sh
#
# Add/Update the header keywords of our sample products corresponding to
# the "repro-4", by running the `r4_header_update' tool.
# These new keywords are required by new version (>= v4.6) tools such as
# `dmcoords', `mkacisrmf', `mkwarf', etc.
#
# Reference:
# [1] ahelp - r4_header_update
#     http://cxc.harvard.edu/ciao/ahelp/r4_header_update.html
#
# Aaron LI
# Created: 2016-03-03
# Updated: 2016-03-03
#


if [ $# -lt 1 ]; then
    echo "Usage:"
    echo "    `basename $0` <repro_dir_list>"
    echo "    `basename $0` <repro_dir1> ..."
    exit 1
fi

# repro dir's
if [ -f "$1" ]; then
    REPRO_LIST="$1"
elif [ -d "$1" ]; then
    REPRO_LIST="_tmp_repro_$$.list"
    [ -f "${REPRO_LIST}" ] && mv -f "${REPRO_LIST}" "${REPRO_LIST}_bak"
    while [ ! -z "$1" ]; do
        echo "$1" >> "${REPRO_LIST}"
        shift
    done
else
    echo "ERROR: invalid arguments: $1"
    exit 2
fi

INIT_DIR=`pwd -P`
cat "${REPRO_LIST}" | while read repro_dir; do
    cd "${repro_dir}"
    echo "********* `pwd` *********"
    PBKFILE=`ls acisf*_pbk0.fits`
    ASOLLIS=`ls acisf*_asol1.lis`
    # fix asol.lis if necessary
    if grep -q '/' "${ASOLLIS}"; then
        sed -i'_bak' -e 's#^/.*/##' "${ASOLLIS}"
    fi
    for f in `find . -type f \( -iname 'acisf*_repro_evt2.fits' -o -iname 'evt2*.fits' -o -iname 'img*.fits' \)`; do
        echo "* ${f}"
        r4_header_update infile="${f}" pbkfile="${PBKFILE}" asolfile="@${ASOLLIS}"
    done
    cd "${INIT_DIR}"
done

