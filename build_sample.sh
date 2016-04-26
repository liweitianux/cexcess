#!/bin/sh
#
# Build and/or update the sample of studying central X-ray brightness
# excess, by only copying/updating the required products from our complete
# mass sample.
#
# Aaron LI
# Created: 2016-03-02
# Updated: 2016-04-15
#
# ChangeLog:
# 2016-04-15:
#   * Add support for zzh's sample
#   * Convert "_" to "-" in NAME
#

analyze_path() {
    # extract `obs_id' and `source name' from path
    #
    # Weitian LI <liweitianux@gmail.com>
    # 2013/02/04
    #
    # input:
    #     path that include `oi' and `source name'
    #     e.g.:
    #
    # output:
    #

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


if [ $# -lt 2 ]; then
    echo "Usage:"
    echo "    `basename $0` <dest_root_dir> <repro_dir_list>"
    echo "    `basename $0` <dest_root_dir> <repro_dir1> ..."
    exit 1
fi

DEST_ROOT="$1"
[ ! -d  "${DEST_ROOT}" ] && mkdir -v "${DEST_ROOT}"

# repro dir's
if [ -f "$2" ]; then
    REPRO_LIST="$2"
elif [ -d "$2" ]; then
    REPRO_LIST="_tmp_repro_$$.list"
    [ -f "${REPRO_LIST}" ] && mv -f "${REPRO_LIST}" "${REPRO_LIST}_bak"
    while [ ! -z "$2" ]; do
        echo "$2" >> "${REPRO_LIST}"
        shift
    done
else
    echo "ERROR: invalid arguments: $2"
    exit 2
fi

INIT_DIR=`pwd -P`
cat "${REPRO_LIST}" | while read repro_dir; do
    OI=`analyze_path "${repro_dir}" | grep '^oi:' | awk '{ print $2 }'`
    NAME=`analyze_path "${repro_dir}" | grep '^name:' | awk '{ print $2 }'`
    OWNER=`analyze_path "${repro_dir}" | grep '^owner:' | awk '{ print $2 }'`
    # avoid "_" in name
    NAME=`echo "${NAME}" | tr '_' '-'`
    echo "********* ${NAME}_oi${OI} *********"
    # create directories
    DEST="${DEST_ROOT}/${NAME}_oi${OI}"
    #[ ! -d "${DEST}/repro" ] && mkdir -pv ${DEST}/repro
    [ ! -d "${DEST}/repro" ] && continue  # Only update sample
    cd "${DEST}/repro"
    # simply copy the whole sub-directories
    cp -av ${repro_dir}/acisf*.fits .
    cp -av ${repro_dir}/acisf*.lis ${repro_dir}/pcadf*.fits .
    cp -av ${repro_dir}/evt ${repro_dir}/bkg ${repro_dir}/img .
    cp -av ${repro_dir}/*_INFO.json .
    if [ "${OWNER}" = "zzh" ]; then
        cp -av ${repro_dir}/spc/profile/rspec.reg ./img/
        cp -av ${repro_dir}/../evt2/spc/profile/rspec.reg ./img/rspec2.reg
        cp -av ${repro_dir}/../evt2/img/sbprofile.reg ./img/sbprofile2.reg
        cd ./img
        [ ! -f "rspec.reg" ] && ln -sv rspec2.reg rspec.reg
        [ ! -f "sbprofile.reg" ] && ln -sv sbprofile2.reg sbprofile.reg
    fi
    # apply clean up
    find . \( -iname '*_bak' -o -iname '_tmp*' \) -delete
    find . -type l -iname '*.dat' -delete
    cd "${INIT_DIR}"
done

