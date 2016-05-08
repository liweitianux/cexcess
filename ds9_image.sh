#!/bin/sh
#
# Use DS9 to visualize the ACIS image, with useful/handy
# arguments/options passed.
# Also touch the output image filename for easier save.
#
# Aaron LI
# Created: 2016-04-16
# Updated: 2016-05-08
#

# Default parameters
FILE_PATTERN="img_c*_fill.fits"
REG_FILE=${REG_FILE:-r500.reg}
REG_FORMAT=${REG_FORMAT:-ds9}
REG_FOV=${REG_FOV:-skyfov.fits}
SMOOTH_RADIUS=${SMOOTH_RADIUS:-4}


case "$1" in
    -[hH]*)
        echo "Usage:"
        echo "    `basename $0` <img_dir1> ..."
        exit 1
        ;;
esac


analyze_path() {
    # extract `obs id' and `source name' from path
    echo "$@" | awk '
    # main part
    {
        if (NF==1) {
            ## oi & name
            input=($1 "/")
            if (input ~ /_oi/) {
                ## PATTERN: .../$name_oi$oi/...
                idx_oi = match(input, /oi[0-9]+/) + 2;    # "2" skip the "oi"
                len_oi = RLENGTH - 2;
                oi = substr(input, idx_oi, len_oi);
                idx_name = match(input, /\/[a-zA-Z0-9.+-]+_oi/) + 1;
                len_name = RLENGTH - 4;
                name = substr(input, idx_name, len_name);
                owner = "lwt";
            }
            else {
                ## PATTERN: .../$name/$oi/...
                idx_oi = match(input, /\/[0-9]+\//) + 1;
                len_oi = RLENGTH - 2;
                oi = substr(input, idx_oi, len_oi);
                idx_name1 = match(input, /\/[a-zA-Z0-9_.+-]+\/[0-9]+\//);
                len_name1 = RLENGTH;
                name1 = substr(input, idx_name1, len_name1);
                idx_name = match(name1, /\/[a-zA-Z0-9_.+-]+\//) + 1;
                len_name = RLENGTH - 2;
                name = substr(name1, idx_name, len_name);
                owner = "zzh";
            }
            ## output
            printf("input: %s\n", input)
            printf("oi: %s\nname: %s\nowner: %s\n", oi, name, owner)
        }
        else {
            printf("*** WARNING: invalid input: %s\n", $0)
        }
    }
    # END { }
    '
}


INIT_DIR=`pwd -P`
while [ ! -z "$1" ]; do
    imgdir="$1"
    shift
    cd ${imgdir}
    echo "====== ${PWD} ======"
    NAME=`analyze_path "${PWD}" | grep '^name:' | awk '{ print $2 }'`
    FILE=`\ls ${FILE_PATTERN}`
    [ -z "${FILE}" ] && continue
    echo "FITS file: ${FILE}"
    #
    PNG_SUFFIX=`echo "${FILE%.fits}" | sed -e 's/_c7//' -e 's/_c0-3//'`
    PNG_FILE="${NAME}_${PNG_SUFFIX}.png"
    [ ! -f "${PNG_FILE}" ] && touch ${PNG_FILE}
    echo "PNG file: ${PNG_FILE}"
    # FoV region
    if [ -f "${REG_FOV}" ]; then
        DS9_REG_FOV="-regions color green -regions ${REG_FOV}"
    fi
    #
    ds9 ${FILE} \
        -bin factor 1 \
        -smooth radius ${SMOOTH_RADIUS} -smooth yes \
        -scale asinh \
        -cmap sls \
        -geometry 1288x1026-0+0 \
        -regions format ${REG_FORMAT} \
        -regions color white -regions ${REG_FILE} ${DS9_REG_FOV}
    #
    cd ${INIT_DIR}
done
