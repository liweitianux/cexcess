#!/bin/sh
#
# Create the new `sbpfit' subdirectory, and prepare the files for fitting
# the surface brightness profile.
#
# Aaron LI
# Created: 2016-03-28
#

prepare() {
    img_dir="$1"
    info="$2"
    ln -sv ${img_dir}/sbprofile.* ${img_dir}/sbprofile_rmid.fits .
    ln -sv ${img_dir}/evt2_c*_clean.fits .
    # sbpfit config
    cp ${SBPFIT_SBETA_CONF} ${SBPFIT_DBETA_CONF} .
    date=`date --iso-8601=seconds`
    sed -i'' -e "s#<DATE>#${date}#" `basename ${SBPFIT_SBETA_CONF}`
    sed -i'' -e "s#<DATE>#${date}#" `basename ${SBPFIT_DBETA_CONF}`
    if [ -n "${info}" ]; then
        name=`grep 'Source Name' ${info} | awk -F'"' '{ print $4 }'`
        obsid=`grep 'Obs. ID' ${info} | awk -F':' '{ print $2 }' | tr -d ' ,'`
        echo "Name: ${name}; ObsID: ${obsid}"
        # sbeta
        sed -i'' -e "s#<NAME>#${name}#"   `basename ${SBPFIT_SBETA_CONF}`
        sed -i'' -e "s#<OBSID>#${obsid}#" `basename ${SBPFIT_SBETA_CONF}`
        # dbeta
        sed -i'' -e "s#<NAME>#${name}#"   `basename ${SBPFIT_DBETA_CONF}`
        sed -i'' -e "s#<OBSID>#${obsid}#" `basename ${SBPFIT_DBETA_CONF}`
    fi
}


if [ $# -ne 2 ]; then
    echo "Usage:"
    echo "    `basename $0` <config_dir> <repro_list>"
    exit 1
fi

CUR_DIR=`pwd -P`
CONFIG_DIR=`realpath $1`
SBPFIT_SBETA_CONF="${CONFIG_DIR}/sbpfit_sbeta.conf"
SBPFIT_DBETA_CONF="${CONFIG_DIR}/sbpfit_dbeta.conf"

cat $2 | while read repro; do
    echo "*** ${repro} ***"
    cd ${CUR_DIR}
    cd ${repro}
    REPRO_DIR=`pwd -P`
    SBPFIT_DIR="${REPRO_DIR}/sbpfit"
    [ -d "${SBPFIT_DIR}" ] && mv -fv ${SBPFIT_DIR} ${SBPFIT_DIR}_bak
    mkdir ${SBPFIT_DIR} && cd ${SBPFIT_DIR}
    if [ -d "../img/ne" ] && [ -d "../img/sw" ]; then
        echo "NOTE: exists 'ne' and 'sw' two parts"
        mkdir ne && cd ne
        INFO=`realpath ${REPRO_DIR}/*ne_INFO.json 2>/dev/null`
        prepare ../../img/ne ${INFO}
        cd ${SBPFIT_DIR}
        mkdir sw && cd sw
        INFO=`realpath ${REPRO_DIR}/*sw_INFO.json 2>/dev/null`
        prepare ../../img/sw ${INFO}
    else
        INFO=`realpath ${REPRO_DIR}/*_INFO.json 2>/dev/null`
        prepare ../img ${INFO}
    fi
done

