#!/bin/sh
#
# Build and/or update the sample of studying central X-ray brightness
# excess, by only copying/updating the required products from our complete
# mass sample.
#
# Aaron LI
# Created: 2016-03-02
# Updated: 2016-06-06
#
# Change log:
# 2016-06-06:
#   * Add argument "OWNER"
#   * Split out function "copy_mass_files()"
#   * Copy the `mass` related files
# 2016-04-26:
#   * Add flag "UPDATE_MODE"
#   * Choose the correct `img` directory for ZZH's data
#   * Simplify the argument processing
# 2016-04-15:
#   * Replace "_" with "-" in NAME
#   * Add support for zzh's sample
#

analyze_path() {
    # extract `obs_id' and `source name' from path
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


# copy mass-related files
copy_mass_files() {
    _mass_dir="$1"
    _target_dir="$2"
    _owner="$3"
    ( \
        cd ${repro_dir}/mass; \
        if [ "${_owner}" = "lwt" ]; then \
            cp -av fitting_*.conf global.cfg wang2012_param.txt \
               flux_sbp.txt radius_sbp.txt sbprofile.txt \
               [a-z]*_center.* summary_*.* final_result.txt \
               ${target_dir}/; \
            cp -Lv tcl_temp_profile.txt ${target_dir}/; \
        else \
            cp -av global.cfg sbp.cfg param_wang.txt \
               [a-z]*_center.* summary_*.* \
               final_result.txt ${target_dir}/; \
            cp -Lv Tprofile.dat radius.dat sbp2.dat sbp3.dat ${target_dir}/; \
        fi; \
    )
    unset _mass_dir _target_dir _owner
}


case "$1" in
    -[hH]*)
        echo "Usage:"
        echo "    `basename $0` <lwt|zzh> <dest_root_dir> <repro_dir1> ..."
        exit 1
        ;;
esac

UPDATE_MODE="YES"

OWNER="$1"
DEST_ROOT="$2"
if [ -n "${UPDATE_MODE}" -a ! -d "${DEST_ROOT}" ]; then
    echo "UPDATE MODE: '${DEST_ROOT}' dose not exist"
    exit 2
fi
[ ! -d  "${DEST_ROOT}" ] && mkdir -v "${DEST_ROOT}"

INIT_DIR=`pwd -P`
while [ ! -z "$3" ]; do
    repro_dir="$3"
    shift
    OI=`analyze_path "${repro_dir}" | grep '^oi:' | awk '{ print $2 }'`
    NAME=`analyze_path "${repro_dir}" | grep '^name:' | awk '{ print $2 }'`
    # avoid "_" in name
    NAME=`echo "${NAME}" | tr '_' '-'`
    echo "********* ${NAME}_oi${OI} *********"
    # create directories
    DEST="${DEST_ROOT}/${NAME}_oi${OI}"
    if [ -n "${UPDATE_MODE}" -a ! -d "${DEST}/repro" ]; then
        # skip if dest repro directory does not exist
        echo "Skipped!"
        continue
    fi
    [ ! -d "${DEST}/repro" ] && mkdir -pv ${DEST}/repro
    cd "${DEST}/repro"
    # copy the re-processed evt2 files
    cp -av ${repro_dir}/acisf*.fits .
    cp -av ${repro_dir}/acisf*.lis ${repro_dir}/pcadf*.fits .
    # copy the INFO.json file if present
    cp -av ${repro_dir}/*_INFO.json .
    # copy the cleaned evt2 files and background files
    cp -av ${repro_dir}/evt .
    cp -av ${repro_dir}/bkg .
    # copy the `img` directory
    if [ "${OWNER}" = "lwt" ]; then
        cp -av ${repro_dir}/img .
    elif [ "${OWNER}" = "zzh" ]; then
        img_dir="${repro_dir}/../evt2/img"
        if [ -f "${img_dir}/sbprofile.txt" ]; then
            cp -av ${img_dir} .
            cp -av ${repro_dir}/../evt2/spc/profile/rspec.reg ./img/
        else
            echo "WARNING: '${img_dir}/sbprofile.txt' does not exists"
            cp -av ${repro_dir}/img .
        fi
    else
        echo "ERROR: unknown source owner: ${OWNER}"
    fi
    # copy the `mass` related files
    mkdir ./mass
    cd ./mass
    target_dir=`pwd -P`
    if [ "${OWNER}" = "lwt" ]; then
        copy_mass_files ${repro_dir}/mass ${target_dir} ${OWNER}
    elif [ "${OWNER}" = "zzh" ]; then
        if [ -d "${mass_dir}/../evt2/mass" ]; then
            copy_mass_files ${repro_dir}/../evt2/mass ${target_dir} ${OWNER}
        elif [ -d "${repro_dir}/mass" ]; then
            copy_mass_files ${repro_dir}/mass ${target_dir} ${OWNER}
        else
            echo "ERROR: no mass directory found!"
        fi
    else
        echo "ERROR: unknown source owner: ${OWNER}"
    fi
    cd ${target_dir}
    # clean up the copied products
    find . \( -iname '*_bak' -o -iname '_tmp*' \) -delete
    find . -type l -iname '*.dat' -delete
    cd ${INIT_DIR}
done
