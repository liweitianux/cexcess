#
# Aaron LI
# Created: 2016-05-04
# Updated: 2016-05-04
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


def get_r500(info):
    """
    Get the R500 value (in unit pixel and kpc), as well as its errors.

    Arguments:
      * info: filename of the INFO.json, or the info dictionary

    Return:
      a dictionary contains the necessary results
    """
    info = load_info(info)

    if "IN_SAMPLE" in info.keys():
        # ZZH's INFO
        r500_kpc = float(info["R500"])
        r500EL_kpc = float(info["R500_err_lower"])
        r500EU_kpc = float(info["R500_err_upper"])
    else:
        # LWT's INFO
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
