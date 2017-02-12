#!/bin/sh
#
# Copyright (c) 2016 Aaron LI
# MIT license
#
# Make the corresponding stowed bacground image for the observation.
#
# Based on `ciao_blanksky.sh' (v5.0; 2015-06-02)
#
# Created: 2016-04-20
#
# Changelog:
# 2016-04-28:
#   * Use the existing only one bgstow if no proper bgstow found
#   * Update keyword "DETNAM" after merge
#

## error code {{{
ERR_USG=1
ERR_DIR=11
ERR_EVT=12
ERR_BKG=13
ERR_REG=14
ERR_ASOL=21
ERR_BPIX=22
ERR_PBK=23
ERR_MSK=24
ERR_BKGTY=31
ERR_SPEC=32
## error code }}}

## usage, help {{{
case "$1" in
    -[hH]*|--[hH]*)
        printf "usage:\n"
        printf "    `basename $0` evt=<evt_fits> basedir=<base_dir> imgdir=<img_dir> erange=<elow:ehigh>\n"
        exit ${ERR_USG}
        ;;
esac
## usage, help }}}

## default parameters {{{
# default `event file' which used to match `blanksky' files
#DFT_EVT="_NOT_EXIST_"
DFT_EVT="`\ls evt2*_clean.fits`"
# default dir which contains `asols, asol.lis, ...' files
# DFT_BASEDIR="_NOT_EXIST_"
DFT_BASEDIR=".."
# default img directory
DFT_IMGDIR="../img"
# default energy range for image creation
DFT_ERANGE="700:7000"

## howto find files in `basedir'
# default `asol.lis pattern'
DFT_ASOLIS_PAT="acis*asol?.lis"
## default parameters }}}

## functions {{{
# process commandline arguments
# cmdline arg format: `KEY=VALUE'
getopt_keyval() {
    until [ -z "$1" ]
    do
        key=${1%%=*}                    # extract key
        val=${1#*=}                     # extract value
        keyval="${key}=\"${val}\""
        echo "## getopt: eval '${keyval}'"
        eval ${keyval}
        shift                           # shift, process next one
    done
}

# lookup the corresponding stowed background according to the blanksky
lookup_bgstow() {
    _dir=`dirname $1`
    _blanksky=`basename $1`
    _det=`echo "${_blanksky}" | sed 's/\(acis[0-7]\).*$/\1/'`
    _year=`echo "${_blanksky}" | sed 's/.*\(D[0-9]\{4\}\).*$/\1/'`
    _cti=`echo "${_blanksky}" | sed 's/.*\(cti\).*$/\1/'`
    if [ "${_year}" = "D1999" ]; then
        _year="D2000"
    fi
    if [ "${_cti}" = "cti" ]; then
        _bgstow=`\ls ${_dir}/${_det}${_year}-??-??bgstow*.fits | grep 'bgstow_cti'`
    else
        _bgstow=`\ls ${_dir}/${_det}${_year}-??-??bgstow*.fits | grep -v 'bgstow_cti'`
    fi
    if [ -z "${_bgstow}" ]; then
        if [ `\ls ${_dir}/${_det}${_year}-??-??bgstow*.fits | wc -l` -eq 1 ]; then
            _bgstow=`\ls ${_dir}/${_det}${_year}-??-??bgstow*.fits`
            echo "WARNING: no proper bgstow; but use the only one" > /dev/stderr
        else
            echo "ERROR: cannot found bgstow for blanksky: ${_dir}/${_blanksky}" >/dev/stderr
            exit 100
        fi
    fi
    echo "${_bgstow}"
    unset _dir _blanksky _det _year _bgstow
}

# reprocess background evt with matched gainfile
background_regain() {
    _outfile="$1"
    _infile="${1%.fits}_ungain.fits"
    mv ${_outfile} ${_infile}
    punlearn acis_process_events
    acis_process_events infile="${_infile}" \
        outfile="${_outfile}" \
        acaofffile=NONE stop="none" doevtgrade=no \
        apply_cti=yes apply_tgain=no \
        calculate_pi=yes pix_adj=NONE \
        gainfile="$2" \
        eventdef="{s:ccd_id,s:node_id,i:expno,s:chip,s:tdet,f:det,f:sky,s:phas,l:pha,l:pha_ro,f:energy,l:pi,s:fltgrade,s:grade,x:status}" \
        clobber=yes
    rm -fv ${_infile}
    unset _infile _outfile
}
## functions end }}}

## parameters {{{
# process cmdline args using `getopt_keyval'
getopt_keyval "$@"

# check given parameters
# check evt file
if [ -r "${evt}" ]; then
    EVT=${evt}
elif [ -r "${DFT_EVT}" ]; then
    EVT=${DFT_EVT}
else
    read -p "evt2 file: " EVT
    if ! [ -r "${EVT}" ]; then
        printf "ERROR: cannot access given \`${EVT}' evt file\n"
        exit ${ERR_EVT}
    fi
fi
printf "## use evt file: \`${EVT}'\n"
# check given dir
if [ -d "${basedir}" ]; then
    BASEDIR=${basedir}
elif [ -d "${DFT_BASEDIR}" ]; then
    BASEDIR=${DFT_BASEDIR}
else
    read -p "basedir (contains asol files): " BASEDIR
    if [ ! -d ${BASEDIR} ]; then
        printf "ERROR: given \`${BASEDIR}' NOT a directory\n"
        exit ${ERR_DIR}
    fi
fi
# remove the trailing '/'
BASEDIR=`echo ${BASEDIR} | sed 's/\/*$//'`
printf "## use basedir: \`${BASEDIR}'\n"
# check img dir
if [ -n "${imgdir}" ]; then
    IMGDIR="${imgdir}"
else
    IMGDIR="${DFT_IMGDIR}"
fi
printf "## use imgdir: \`${IMGDIR}'\n"
# check energy range
if [ -n "${erange}" ]; then
    ERANGE="${erange}"
else
    ERANGE="${DFT_ERANGE}"
fi
printf "## use energy range: \`${ERANGE}'\n"
## parameters }}}

## check files in `basedir' {{{
# check asol files
ASOLIS=`\ls -1 ${BASEDIR}/${DFT_ASOLIS_PAT} | head -n 1`
if [ -z ${ASOLIS} ]; then
    printf "ERROR: cannot find \"${DFT_ASOLIS_PAT}\" in dir \`${BASEDIR}'\n"
    exit ${ERR_ASOL}
fi
printf "## use asolis: \`${ASOLIS}'\n"
## check files }}}

## prepare parameter files (pfiles) {{{
CIAO_TOOLS="acis_bkgrnd_lookup dmmerge dmcopy dmmakepar dmreadpar reproject_events acis_process_events"

# Copy necessary pfiles for localized usage
for tool in ${CIAO_TOOLS}; do
    pfile=`paccess ${tool}`
    [ -n "${pfile}" ] && punlearn ${tool} && cp -Lvf ${pfile} .
done

# Modify environment variable 'PFILES' to use local pfiles first
export PFILES="./:${PFILES}"
## pfiles }}}

#### main start {{{
# decompress the asol file if necessary
ASOL_XZ=`\ls ${BASEDIR}/pcadf*asol?.fits.xz 2>/dev/null`
if [ -n "${ASOL_XZ}" ]; then
    printf "decompressing the asol file ...\n"
    unxz ${ASOL_XZ}
fi

printf "look up coresponding background file ...\n"
punlearn acis_bkgrnd_lookup
BKG_LKP="`acis_bkgrnd_lookup ${EVT}`"
AS_NUM=`echo ${BKG_LKP} | tr ' ' '\n' | \grep 'acis7sD' | wc -l`
AI_NUM=`echo ${BKG_LKP} | tr ' ' '\n' | \grep 'acis[0123]iD' | wc -l`
## determine detector type: ACIS-S / ACIS-I {{{
if [ ${AS_NUM} -eq 1 ]; then
    printf "## ACIS-S, chip: 7\n"
    CCD="7"
    DETNAM="ACIS-7"
    BKG_ROOT="stowbkg_c7"
    STOWBKG=`lookup_bgstow ${BKG_LKP}`
    cp -v ${STOWBKG} ${BKG_ROOT}_orig.fits
elif [ ${AI_NUM} -eq 4 ]; then
    printf "## ACIS-I, chip: 0-3\n"
    CCD="0:3"
    DETNAM="ACIS-0123"
    BKG_ROOT="stowbkg_c0-3"
    AI_FILES=""
    for bkg in ${BKG_LKP}; do
        STOWBKG=`lookup_bgstow ${bkg}`
        cp -v ${STOWBKG} .
        AI_FILES="${AI_FILES},`basename ${STOWBKG}`"
    done
    AI_FILES=${AI_FILES#,}      # remove the first ','
    printf "## ACIS-I background files to merge:\n"
    printf "##   \`${AI_FILES}'\n"
    printf "\`dmmerge' to merge the above blanksky files ...\n"
    # merge 4 chips blanksky evt files
    punlearn dmmerge
    dmmerge "${AI_FILES}" ${BKG_ROOT}_orig.fits clobber=yes
    # update DETNAM
    punlearn dmhedit
    dmhedit infile="${BKG_ROOT}_orig.fits" filelist=none operation=add \
        key=DETNAM value="${DETNAM}"
    rm -fv `echo ${AI_FILES} | tr ',' ' '`      # remove original files
else
    printf "## ERROR: UNKNOW blanksky files:\n"
    printf "##   ${BKG_ORIG}\n"
    exit ${ERR_BKG}
fi
## determine ACIS type }}}

## check 'DATMODE' {{{
## filter blanksky files (status=0) for `VFAINT' observations
DATA_MODE="`dmkeypar ${EVT} DATAMODE echo=yes`"
printf "## DATAMODE: ${DATA_MODE}\n"
if [ "${DATA_MODE}" = "VFAINT" ]; then
    mv -fv ${BKG_ROOT}_orig.fits ${BKG_ROOT}_tmp.fits
    printf "apply \`status=0' to filter blanksky file ...\n"
    punlearn dmcopy
    dmcopy "${BKG_ROOT}_tmp.fits[status=0]" ${BKG_ROOT}_orig.fits clobber=yes
    rm -fv ${BKG_ROOT}_tmp.fits
fi
## DATAMODE, status=0 }}}

## check `GAINFILE' of blanksky and evt2 file {{{
## if NOT match, then reprocess blanksky
GAINFILE_EVT="`dmkeypar ${EVT} GAINFILE echo=yes`"
GAINFILE_BG="`dmkeypar "${BKG_ROOT}_orig.fits" GAINFILE echo=yes`"
if ! [ "${GAINFILE_EVT}" = "${GAINFILE_BG}" ]; then
    printf "WARNING: GAINFILE NOT match.\n"
    printf "event: ${GAINFILE_EVT}\n"
    printf "blank: ${GAINFILE_BG}\n"
    printf "reprocess blanksky with evt gainfile ...\n"
    # reprocess blanksky using matched evt GAINFILE
    GAINFILE="$CALDB/data/chandra/acis/det_gain/`basename ${GAINFILE_EVT}`"
    printf "GAINFILE: ${GAINFILE}\n"
    background_regain "${BKG_ROOT}_orig.fits" ${GAINFILE}
fi
## check & match GAINFILE }}}

printf "add the PNT header keywords ... "
EVT_HEADER="_evt_header.par"
EVT_PNT="_evt_pnt.par"
punlearn dmmakepar
dmmakepar ${EVT} ${EVT_HEADER} clobber=yes
\grep -i '_pnt' ${EVT_HEADER} > ${EVT_PNT}
punlearn dmreadpar
dmreadpar ${EVT_PNT} "${BKG_ROOT}_orig.fits[EVENTS]" clobber=yes
printf "DONE\n"

printf "reproject the background ...\n"
punlearn reproject_events
reproject_events infile=${BKG_ROOT}_orig.fits \
    outfile=${BKG_ROOT}.fits match=${EVT} \
    aspect="@${ASOLIS}" random=0 clobber=yes

# Make image of the stowed background
if [ -f "${IMGDIR}/skyfov.fits" ]; then
    ln -svf ${IMGDIR}/skyfov.fits .
    BKG_IMG="img_${BKG_ROOT}_e`echo ${ERANGE} | tr ':' '-'`.fits"
    printf "creating the background image ...\n"
    punlearn dmcopy
    dmcopy infile="${BKG_ROOT}.fits[energy=${ERANGE}][sky=region(skyfov.fits[ccd_id=${CCD}])][bin sky=1]" outfile=${BKG_IMG} clobber=yes
else
    printf "ERROR: '${IMGDIR}/skyfov.fits' not exists\n"
fi
## main end }}}

# clean
printf "\nclean ...\n"
rm -fv ${BKG_ROOT}_orig.fits ${EVT_PNT}

printf "\nFINISHED\n"

# vim: set ts=8 sw=4 tw=0 fenc=utf-8 ft=sh #
