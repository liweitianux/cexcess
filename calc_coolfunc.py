#!/usr/bin/env python3
#
# Based on 'coolfunc_calc.sh' from 'chandra-acis-analysis/mass_profile'.
#
# Aaron LI
# Created: 2016-06-19
# Updated: 2016-07-04
#
# Change logs:
# 2016-07-04:
#   * Use "AstroParams"
#   * Update to use 3-column temperature profile
#   * Update documentation
#   * Update sample configuration file
# 2016-06-27:
#   * Minor style fixes
#   * Change 'tprofile' to 't_profile'
#

"""
Calculate the *cooling function* profile with respect to the
given *temperature profile* and the average abundance, redshift,
and column density nH, using the XSPEC model 'wabs*apec'.

Emission measure:
    EM = \int n_e n_H dV ~= (n_e^2 / ratio_eH) V  [ cm^-3 ]
where 'ratio_eH' is the ratio of electron density to proton density (n_H).

APEC normalization returned by XSPEC is simply the *emission measure* of
the gas scaled by the distance:
    eta = (\int n_e n_H dV) / (4 pi (D_A (1+z))^2)

The flux calculated with the XSPEC `flux` command has dimension:
    Flux: [ photon s^-1 cm^-2 ] or [ erg s^-2 cm^-2 ]

If we let EM=1 and then set the APEC's normalization,
the cooling function is therefore derived by calculating the flux
using the XSPEC `flux` command, and the cooling function has
dimension of [ FLUX / EM ].

See also the documentation of `deproject_sbp.py` for more details.


Sample configuration file:
------------------------------------------------------------
## Configuration file for `calc_coolfunc.py`
## 2016-07-04

# temperature profile fitted & extrapolated by model: [r, T]
t_profile = t_profile.txt

# average abundance (unit: solar)
abundance = 0.5

# abundance table (default: grsa)
abund_table = grsa

# redshift of the object
redshift = <REDSHIFT>

# H column density (unit: 10^22 cm^-2)
nh = <NH>

# energy range within which to calculate the cooling function (unit: keV)
energy_low = 0.7
energy_high = 7.0

# output file of the XSPEC script for cooling function calculation
xspec_script = coolfunc.xcm

# output file of the cooling function profile: [r, CF]
coolfunc = coolfunc_profile.txt
------------------------------------------------------------
"""

import argparse
import subprocess
import sys
import os
from datetime import datetime

import numpy as np
import astropy.units as au
from astropy.cosmology import FlatLambdaCDM
from configobj import ConfigObj

from astro_params import AstroParams


def gen_xspec_script(outfile, data):
    """
    Generate the XSPEC script for cooling function profile calculation.

    Arguments:
    * outfile: output file to save the XSPEC script
    * data: dictionary used to format the template XSPEC script
    """
    xspec_script = """
# Calculate the cooling function profile w.r.t the temperature profile.
#
# Generated by: %(prog_name)s
# Date: %(cur_date)s

# debug (off)
chatter 0

set xs_return_results 1
set xs_echo_script 0
# set tcl_precision 12

query yes
abund %(abund_table)s
dummyrsp 0.01 100.0 4096 linear
# use model 'wabs*apec'
model wabs*apec & %(nh)s & 1.0 & %(abundance)s & %(redshift)s & %(apec_norm)s & /*

# input and output files
set tpro_fn "%(t_profile)s"
set cf_fn "%(coolfunc)s"
if { [ file exists $cf_fn ] } {
    exec rm -fv $cf_fn
}

# open files
set tpro_fd [ open $tpro_fn r ]
set cf_fd [ open $cf_fn w ]

# output file header
puts $cf_fd "# radius    flux(%(energy_low)s-%(energy_high)s)"

# read data from temperature profile line by line
while { [ gets $tpro_fd line ] != -1 } {
    if {[ regexp -- {^\s*#} $line ] == 1} {
        # ignore comment line
        continue
    }
    scan $line "%%f %%f %%f" radius radius_err temperature
    #puts "radius: $radius, temperature: $temperature
    # set temperature value
    newpar 2 $temperature
    flux %(energy_low)s %(energy_high)s
    tclout flux 1
    scan $xspec_tclout "%%f %%f %%f %%f" _ _ _ cf_data
    #puts "cf_data: $cf_data"
    puts $cf_fd "$radius    $cf_data"
}

# close & exit
close $tpro_fd
close $cf_fd
tclexit
""" % data
    open(outfile, "w").write(xspec_script)


def calc_apec_norm(z):
    """
    Calculate the normalization of the APEC model.

    Reference:
    https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/XSmodelApec.html
    """
    cosmo = FlatLambdaCDM(H0=AstroParams.H0, Om0=AstroParams.OmegaM0)
    D_A = cosmo.angular_diameter_distance(z).to(au.cm).value
    norm = 1.0e-14 / (4*np.pi * (D_A * (1+z))**2)
    return norm


def calc_coolfunc(xspec_script, verbose=True):
    if verbose:
        print("Invoke XSPEC to calculate cooling function profile ...")
    subprocess.run(args=["xspec", "-", xspec_script],
                   stdout=subprocess.DEVNULL)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate the cooling function profile " +
                    "w.r.t the temperature profile")
    parser.add_argument("config", nargs="?", default="coolfunc.conf",
                        help="config for cooling function calculation " +
                             "(default: coolfunc.conf)")
    args = parser.parse_args()
    config = ConfigObj(args.config)

    redshift = config.as_float("redshift")
    config_data = {
        "prog_name":    os.path.basename(sys.argv[0]),
        "cur_date":     datetime.now().isoformat(),
        #
        "t_profile":    config["t_profile"],
        "abundance":    config.as_float("abundance"),
        "abund_table":  config.get("abund_table", "grsa"),
        "redshift":     redshift,
        "nh":           config.as_float("nh"),
        "energy_low":   float(config.get("energy_low", 0.7)),
        "energy_high":  float(config.get("energy_high", 0.7)),
        "xspec_script": config["xspec_script"],
        "coolfunc":     config["coolfunc"],
        "apec_norm":    calc_apec_norm(z=redshift),
    }

    gen_xspec_script(outfile=config["xspec_script"], data=config_data)
    calc_coolfunc(xspec_script=config["xspec_script"])


if __name__ == "__main__":
    main()
