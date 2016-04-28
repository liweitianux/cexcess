#!/bin/sh
#
# Generate the exposure map and apply exposure correction
# using CIAO `fluximage`.
#
# NOTE:
# The existing "instmap_weights.txt" is ued for exposure map generation.
#
# Aaron LI
# Created: 2016-04-28
# Updated: 2016-04-28
#

case "$1" in
    -[hH]*)
        echo "Usage:"
        echo "    `basename $0` <evt_e> <img> [ instmap_weights.txt ]"
        exit 1
        ;;
esac

EVT_E="$1"
IMG="$2"
SPEC_WGT="${3:-instmap_weights.txt}"

REPRO_DIR=".."
ASOLIS=`\ls ${REPRO_DIR}/acisf*asol?.lis`
BPIX=`\ls ${REPRO_DIR}/acisf*repro_bpix?.fits`
MSK=`\ls ${REPRO_DIR}/acisf*msk?.fits`

ROOTNAME=`echo "${EVT_E%.fits}" | sed -e 's/^evt2_//'`

## get `xygrid' for image
punlearn get_sky_limits
get_sky_limits image="${IMG}" verbose=0
XYGRID=`pget get_sky_limits xygrid`
echo "xygrid: ${XYGRID}"

punlearn ardlib

echo "invoking fluximage to generate expmap and apply correction ..."
punlearn fluximage
fluximage infile="${EVT_E}" outroot="${ROOTNAME}" \
    binsize=1 bands="${SPEC_WGT}" xygrid="${XYGRID}" \
    asol="@${ASOLIS}" badpixfile="${BPIX}" \
    maskfile="${MSK}" clobber=yes

# make symbolic links
# clipped counts image
ln -svf ${ROOTNAME}*band*thresh.img img_${ROOTNAME}_thresh.fits
# clipped exposure map
ln -svf ${ROOTNAME}*band*thresh.expmap expmap_${ROOTNAME}.fits
# exposure-corrected image
ln -svf ${ROOTNAME}*band*flux.img img_expcorr_${ROOTNAME}.fits

