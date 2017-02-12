#!/bin/sh
#
# Copyright (c) 2016 Aaron LI
# MIT license
#
# Collect the fitting results and orgnize in CSV format.
#
# Created: 2016-03-29
#

if [ $# -ne 1 ]; then
    echo "Usage:"
    echo "    `basename $0` <sbpfit_dir_list>"
fi

SBPFIT_CONF="sbpfit_sbeta.conf"
SBPFIT_DBETA_CONF="sbpfit_dbeta.conf"
CUR=`pwd -P`

# header
printf "Name,ObsID,npoints,"
printf "sbeta_ignore,sbeta_npoints,sbeta_chisq,sbeta_redchisq,s0,s0_E68L,s0_E68U,rc,rc_E68L,rc_E68U,beta,beta_E68L,beta_E68U,bkg,bkg_E68L,bkg_E68U,"
printf "dbeta_ignore,dbeta_npoints,dbeta_chisq,dbeta_redchisq,s01,s01_E68L,s01_E68U,rc1,rc1_E68L,rc1_E68U,beta1,beta1_E68L,beta1_E68U,s02,s02_E68L,s02_E68U,rc2,rc2_E68L,rc2_E68U,beta2,beta2_E68L,beta2_E68U,dbeta_bkg,dbeta_bkg_E68L,dbeta_bkg_E68U\n"

cat $1 | while read sbpfit_dir; do
    cd ${CUR}
    cd ${sbpfit_dir}
    NAME=`grep '^name' ${SBPFIT_CONF} | awk -F'=' '{ print $2 }' | xargs`
    OBSID=`grep '^obsid' ${SBPFIT_CONF} | awk -F'=' '{ print $2 }' | xargs`
    SBP_FILE=`grep '^sbpfile' ${SBPFIT_CONF} | awk -F'=' '{ print $2 }' | xargs`
    RES_FILE=`grep '^outfile' ${SBPFIT_CONF} | awk -F'=' '{ print $2 }' | xargs`
    DBETA_RES_FILE=`grep '^outfile' ${SBPFIT_DBETA_CONF} | awk -F'=' '{ print $2 }' | xargs`
    npoints=`grep -vE '^\s*#' ${SBP_FILE} | wc -l`
    ## single-beta model fitting results
    ignore=`grep '^ignore' ${SBPFIT_CONF} | awk -F'=' '{ print $2 }' | xargs | tr ',' ';'`
    npoints_fit=`grep -E '^\s+# data points' ${RES_FILE} | awk -F'=' '{ print $2 }' | xargs`
    chisq=`grep -E '^\s+chi-square' ${RES_FILE} | awk -F'=' '{ print $2 }' | xargs`
    redchisq=`grep -E '^\s+reduced chi-square' ${RES_FILE} | awk -F'=' '{ print $2 }' | xargs`
    s0=`grep   '^s0:'   ${RES_FILE} | awk '{ print $4 }'`
    rc=`grep   '^rc:'   ${RES_FILE} | awk '{ print $4 }'`
    beta=`grep '^beta:' ${RES_FILE} | awk '{ print $4 }'`
    bkg=`grep  '^bkg:'  ${RES_FILE} | awk '{ print $4 }'`
    # 68% confidence interval
    s0_E68L=`grep   '^s0:'   ${RES_FILE} | awk '{ print $3 }'`
    s0_E68U=`grep   '^s0:'   ${RES_FILE} | awk '{ print $5 }'`
    rc_E68L=`grep   '^rc:'   ${RES_FILE} | awk '{ print $3 }'`
    rc_E68U=`grep   '^rc:'   ${RES_FILE} | awk '{ print $5 }'`
    beta_E68L=`grep '^beta:' ${RES_FILE} | awk '{ print $3 }'`
    beta_E68U=`grep '^beta:' ${RES_FILE} | awk '{ print $5 }'`
    bkg_E68L=`grep  '^bkg:'  ${RES_FILE} | awk '{ print $3 }'`
    bkg_E68U=`grep  '^bkg:'  ${RES_FILE} | awk '{ print $5 }'`
    printf "${NAME},${OBSID},${npoints},"
    printf "${ignore},${npoints_fit},${chisq},${redchisq},${s0},${s0_E68L},${s0_E68U},${rc},${rc_E68L},${rc_E68U},${beta},${beta_E68L},${beta_E68U},${bkg},${bkg_E68L},${bkg_E68U},"
    if [ -r "${DBETA_RES_FILE}" ]; then
        ## double-beta model fitting results
        dbeta_ignore=`grep '^ignore' ${SBPFIT_DBETA_CONF} | awk -F'=' '{ print $2 }' | xargs | tr ',' ';'`
        dbeta_npoints_fit=`grep -E '^\s+# data points' ${RES_FILE} | awk -F'=' '{ print $2 }' | xargs`
        dbeta_chisq=`grep -E '^\s+chi-square' ${DBETA_RES_FILE} | awk -F'=' '{ print $2 }' | xargs`
        dbeta_redchisq=`grep -E '^\s+reduced chi-square' ${DBETA_RES_FILE} | awk -F'=' '{ print $2 }' | xargs`
        s01=`grep       '^s01:'   ${DBETA_RES_FILE} | awk '{ print $4 }'`
        rc1=`grep       '^rc1:'   ${DBETA_RES_FILE} | awk '{ print $4 }'`
        beta1=`grep     '^beta1:' ${DBETA_RES_FILE} | awk '{ print $4 }'`
        s02=`grep       '^s02:'   ${DBETA_RES_FILE} | awk '{ print $4 }'`
        rc2=`grep       '^rc2:'   ${DBETA_RES_FILE} | awk '{ print $4 }'`
        beta2=`grep     '^beta2:' ${DBETA_RES_FILE} | awk '{ print $4 }'`
        dbeta_bkg=`grep '^bkg:'   ${DBETA_RES_FILE} | awk '{ print $4 }'`
        # 68% confidence interval
        s01_E68L=`grep       '^s01:'   ${DBETA_RES_FILE} | awk '{ print $3 }'`
        s01_E68U=`grep       '^s01:'   ${DBETA_RES_FILE} | awk '{ print $5 }'`
        rc1_E68L=`grep       '^rc1:'   ${DBETA_RES_FILE} | awk '{ print $3 }'`
        rc1_E68U=`grep       '^rc1:'   ${DBETA_RES_FILE} | awk '{ print $5 }'`
        beta1_E68L=`grep     '^beta1:' ${DBETA_RES_FILE} | awk '{ print $3 }'`
        beta1_E68U=`grep     '^beta1:' ${DBETA_RES_FILE} | awk '{ print $5 }'`
        s02_E68L=`grep       '^s02:'   ${DBETA_RES_FILE} | awk '{ print $3 }'`
        s02_E68U=`grep       '^s02:'   ${DBETA_RES_FILE} | awk '{ print $5 }'`
        rc2_E68L=`grep       '^rc2:'   ${DBETA_RES_FILE} | awk '{ print $3 }'`
        rc2_E68U=`grep       '^rc2:'   ${DBETA_RES_FILE} | awk '{ print $5 }'`
        beta2_E68L=`grep     '^beta2:' ${DBETA_RES_FILE} | awk '{ print $3 }'`
        beta2_E68U=`grep     '^beta2:' ${DBETA_RES_FILE} | awk '{ print $5 }'`
        dbeta_bkg_E68L=`grep '^bkg:'   ${DBETA_RES_FILE} | awk '{ print $3 }'`
        dbeta_bkg_E68U=`grep '^bkg:'   ${DBETA_RES_FILE} | awk '{ print $5 }'`
        printf "${dbeta_ignore},${dbeta_npoints_fit},${dbeta_chisq},${dbeta_redchisq},${s01},${s01_E68L},${s01_E68U},${rc1},${rc1_E68L},${rc1_E68U},${beta1},${beta1_E68L},${beta1_E68U},${s02},${s02_E68L},${s02_E68U},${rc2},${rc2_E68L},${rc2_E68U},${beta2},${beta2_E68L},${beta2_E68U},${dbeta_bkg},${dbeta_bkg_E68L},${dbeta_bkg_E68U}\n"
    else
        printf ",,,,,,,,,,,,,,,,,,,,,,\n"
    fi
done
