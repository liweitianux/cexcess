#
# Aaron LI
# Created: 2016-05-04
# Updated: 2016-07-12
#
# Change logs:
# 2016-07-12:
#   * Add functions "get_name()" and "get_obsid()"
# 2016-07-11:
#   * Add functions "get_redshift()" and "get_nh()"
#   * Add function "get_owner()"
#

"""
module to process the INFO.json file, contains some handy functions
to extract the needed information.
"""

import json


def load_info(info):
    """
    Load data from the INFO.json if necessary.
    """
    if isinstance(info, str):
        json_str = open(info).read().rstrip().rstrip(",")
        return json.loads(json_str)
    else:
        # assuming that the provided `info` is already a Python dictionary
        return info


def get_owner(info):
    """
    Determine the owner of the info: 'LWT' or 'ZZH'.

    Return:
    * 'LWT'
    * 'ZZH'
    """
    info = load_info(info)
    if "IN_SAMPLE" in info.keys():
        return "ZZH"
    else:
        return "LWT"


def get_name(info):
    """
    Get the source name from the INFO json file.
    """
    info = load_info(info)
    name = info["Source Name"]
    uname = info.get("Unified Name")
    results = {
        "name": name,
        "uname": uname,
    }
    return results


def get_obsid(info):
    """
    Get the Chandra observation ID from the INFO json file.
    """
    info = load_info(info)
    obsid = int(info["Obs. ID"])
    return obsid


def get_r500(info):
    """
    Get the R500 value (in unit pixel and kpc), as well as its errors.

    Arguments:
      * info: filename of the INFO.json, or the info dictionary

    Return:
      a dictionary contains the necessary results
    """
    info = load_info(info)

    if get_owner(info) == "ZZH":
        r500_kpc = float(info["R500"])
        r500EL_kpc = float(info["R500_err_lower"])
        r500EU_kpc = float(info["R500_err_upper"])
    else:
        r500_kpc = float(info["R500 (kpc)"])
        r500EL_kpc = float(info["R500_err_lower (1sigma)"])
        r500EU_kpc = float(info["R500_err_upper (1sigma)"])

    # Convert kpc to Chandra ACIS pixel
    rmax_sbp_pix = float(info["Rmax_SBP (pixel)"])
    rmax_sbp_kpc = float(info["Rmax_SBP (kpc)"])
    kpc_per_pix = rmax_sbp_kpc / rmax_sbp_pix
    r500_pix = r500_kpc / kpc_per_pix
    r500EL_pix = r500EL_kpc / kpc_per_pix
    r500EU_pix = r500EU_kpc / kpc_per_pix

    results = {
            "r500_kpc":    r500_kpc,
            "r500EL_kpc":  r500EL_kpc,
            "r500EU_kpc":  r500EU_kpc,
            "r500_pix":    r500_pix,
            "r500EL_pix":  r500EL_pix,
            "r500EU_pix":  r500EU_pix,
            "kpc_per_pix": kpc_per_pix,
    }
    return results


def get_redshift(info):
    """
    Get the redshift from the info json file.
    """
    info = load_info(info)
    redshift = float(info["redshift"])
    return redshift


def get_nh(info):
    """
    Get the column density (nH) from the info json file.
    """
    info = load_info(info)
    if get_owner(info) == "LWT":
        nh = float(info["nH (10^22 cm^-2)"])
    else:
        nh = float(info["nh"])
    return nh
