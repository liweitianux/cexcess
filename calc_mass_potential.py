#!/usr/bin/env python3
#
# Weitian LI
# Created: 2016-06-24
# Updated: 2016-06-24
#
# Change logs:
# 2016-06-24:
#   * Update method 'gen_radius()'
#

"""
Calculate the (gas and gravitational) mass profile and gravitational
potential profile from the electron number density profile.
The temperature profile is required.

References:
[1] Ettori et al, 2013, Space Science Review, 177, 119-154


Sample configuration file:
------------------------------------------------------------
## Configuration for `calc_mass_potential.py`
## Date: 2016-06-24

# redshift used for pixel to distance conversion
redshift = <REDSHIFT>

# electron density profile
ne_profile = ne_profile.txt

# cooling function profile
cf_profile = coolfunc_profile.txt
# unit of the CF profile radius (default: pixel)
cf_unit = "pixel"

# temperature profile
t_profile = t_profile.txt
# unit of the T profile radius (default: pixel)
t_unit = "pixel"

# number of data points for the output profile calculation
num_dp = 100

# output gas mass profile
m_gas_profile = mass_gas_profile.txt

------------------------------------------------------------
"""

import argparse

import numpy as np
import astropy.units as au
import scipy.interpolate as interpolate
import scipy.integrate as integrate
from configobj import ConfigObj

from astro_params import AstroParams, ChandraPixel
from projection import Projection


class DensityProfile:
    """
    Utilize the 3D (electron number or gas mass) density profile to
    calculate the following quantities:
    * 2D projected surface brightness (requires cooling function profile)
    * gas mass profile (integrated, M_gas(<r))
    * gravitational mass profile (M(<r); requires temperature profile)
    * gravitational potential profile (cut at the largest available radius)

    NOTE:
    * The radii (of density profile and cooling function profile)
      should have unit [ cm ]
    * The density should have unit [ cm^-3 ] or [ g cm^-3 ]
    """
    # allowed density profile types
    DENSITY_TYPES = ["electron", "gas"]
    # input density data: [r, r_err, d]
    r = None
    r_err = None
    d = None
    # electron number density
    ne = None
    # gas mass density
    rho_gas = None
    # cooling function profile
    cf_radius = None
    cf_value = None
    # temperature profile
    t_radius = None
    t_value = None
    # interpolated profiles
    d_interp = None
    cf_interp = None
    t_interp = None
    # generated radial data points for profile calculation
    radius = None
    radius_err = None
    # gas mass profile: M_gas(<r); same length as the above 'radius'
    m_gas = None
    # total (gravitational) mass profile: M_total(<r)
    m_total = None
    # potential profile (cut at the largest available radius)
    potential = None

    def __init__(self, density, density_type="electron"):
        self.load_data(data=density, density_type=density_type)

    def load_data(self, data, density_type="electron"):
        if density_type not in self.DENSITY_TYPES:
            raise ValueError("invalid density_types: %s" % density_type)
        # 3-column density profile: [r, r_err, density]
        self.r = data[:, 0].copy()
        self.r_err = data[:, 1].copy()
        self.d = data[:, 2].copy()
        self.density_type = density_type

    def load_cf_profile(self, data):
        # 2-column cooling function profile: r[cm], cf[flux/EM]
        self.cf_radius = data[:, 0].copy()
        self.cf_value = data[:, 1].copy()

    def load_t_profile(self, data):
        # 2-column temperature profile: r[cm], T[keV]
        self.t_radius = data[:, 0].copy()
        self.t_value = data[:, 1].copy()

    def calc_brightness(self):
        """
        Project the electron number density or gas mass density profile
        to calculate the 2D surface brightness profile.
        """
        if self.cf_radius is None or self.cf_value is None:
            raise ValueError("cooling function profile missing")
        ne = self.calc_electron_density()
        # flux per unit volume
        flux = self.cf_interp(self.r) * ne ** 2 / AstroParams.ratio_ne_np
        # project the 3D flux
        projector = Projection(rout=self.r+self.r_err)
        brightness = projector.project(flux)
        return brightness

    def calc_electron_density(self):
        """
        Calculate the electron number density from the gas mass density
        if necessary.
        """
        if self.density_type == "electron":
            self.ne = self.d.copy()
        elif self.density_type == "gas":
            self.ne = self.d / AstroParams.m_atom / AstroParams.mu_e
        return self.ne

    def calc_gas_density(self):
        """
        Calculate the gas mass density from the electron number density
        if necessary.
        """
        if self.density_type == "electron":
            self.rho_gas = self.d * AstroParams.mu_e * AstroParams.m_atom
        elif self.density_type == "gas":
            self.rho_gas = self.d.copy()
        return self.rho_gas

    def interpolate(self):
        """
        Interpolate the density profile, cooling function profile,
        and temperature profile.

        NOTE:
        * Linear interpolation may be bad because the total mass calculation
          needs to take the derivative of electron density profile and
          temperature profile.  Therefore, smooth spline is a better choice.
        * Allow cooling function profile and temperature profile to be
          extrapolated by filling with the last value if necessary.

        XXX:
        * How to extrapolate the smooth spline if necessary ??
        """
        # density profile
        # insert a data point at radius of zero
        self.d_interp = interpolate.interp1d(
            x=np.concatenate([[0.0], self.r]),
            y=np.concatenate([[self.d[0]], self.d]),
            kind="linear",
            bounds_error=False, fill_value=self.d[-1],
            assume_sorted=True)
        if self.ne is not None:
            self.ne_interp = interpolate.interp1d(
                x=np.concatenate([[0.0], self.r]),
                y=np.concatenate([[self.ne[0]], self.ne]),
                kind="linear",
                bounds_error=False, fill_value=self.ne[-1],
                assume_sorted=True)
        if self.rho_gas is not None:
            self.rho_gas_interp = interpolate.interp1d(
                x=np.concatenate([[0.0], self.r]),
                y=np.concatenate([[self.rho_gas[0]], self.rho_gas]),
                kind="linear",
                bounds_error=False, fill_value=self.rho_gas[-1],
                assume_sorted=True)
        # cooling function profile
        self.cf_interp = interpolate.interp1d(
            x=self.cf_radius, y=self.cf_value, kind="linear",
            bounds_error=False, fill_value=self.cf_value[-1],
            assume_sorted=True)
        # temperature profile
        self.t_interp = interpolate.interp1d(
            x=self.t_radius, y=self.t_value, kind="linear",
            bounds_error=False, fill_value=self.t_value[-1],
            assume_sorted=True)

    def gen_radius(self, num=1000):
        """
        Generate radial points for following mass and potential calculation.

        The generated radial points are logarithmic-evenly distributed.

        NOTE:
        The radii are first generated to determine the inner-most bin width,
        which is used to further divide the original inner-most bin (i.e.,
        radius 0 - r_out[0]), and then the other radii are generated with
        the constraint of given total number of points.
        """
        rout = self.r + self.r_err
        # first pass to determine the inner-most bin width
        rout_tmp = np.logspace(np.log10(rout[0]), np.log10(rout[-1]),
                               num=num, base=10.0)
        bw = rout_tmp[1] - rout_tmp[0]
        # linear-evenly divide the first original bin (0 - rout[0])
        nbin = int(np.ceil(rout[0] / bw))
        rout_new1 = np.linspace(0.0, rout[0], num=nbin, endpoint=False)[1:]
        # second pass to generate the other radii
        rout_new2 = np.logspace(np.log10(rout[0]), np.log10(rout[-1]),
                                num=(num-nbin+1), base=10.0)
        rout_new = np.concatenate([rout_new1, rout_new2])
        rin_new = np.concatenate([[0.0], rout_new[:-1]])
        self.radius = (rout_new + rin_new) / 2.0
        self.radius_err = (rout_new - rin_new) / 2.0

    def calc_mass_gas(self, verbose=False):
        """
        Calculate the gas mass profile, i.e., the mass of the gas within
        each radius.

        Reference: ref.[1], eq.(9)
        """
        def _f_rho_gas(r):
            return self.rho_gas_interp(r) * 4*np.pi * r**2
        #
        m_gas = np.zeros(self.radius.shape)
        if verbose:
            print("Calculating the gas mass profile (#%d): ..." % len(m_gas),
                  end="", flush=True)
        for i, r in enumerate(self.radius):
            if verbose and (i+1) % 10 == 0:
                print("%d..." % (i+1), end="", flush=True)
            # integrated gas mass [ g ]
            m_gas[i] = integrate.quad(_f_rho_gas, a=0.0, b=r,
                                      epsabs=1.0e5, epsrel=1.e-2)[0]
        if verbose:
            print("DONE!", flush=True)
        self.m_gas = m_gas
        return m_gas

    def calc_mass_total(self, verbose=True):
        """
        Calculate the total mass (i.e., gravitational mass) profile,
        under the assumption of hydrostatic equilibrium (HE).

        References: ref.[1], eq.(5,6,7)
        """
        if self.cf_radius is None or self.cf_value is None:
            raise ValueError("cooling function profile required")
        if self.t_radius is None or self.t_value is None:
            raise ValueError("temperature profile required")
        #
        m_total = np.zeros(self.radius.shape)
        if verbose:
            print("Calculating the total mass profile (#%d) ... " %
                  len(m_total))
            # TODO:
        self.m_total = m_total
        return m_total

    def save(self, profile, outfile):
        if profile == "mass_gas":
            data = np.column_stack([self.radius,
                                    self.radius_err,
                                    self.m_gas])
            header = "radius[cm]  radius_err[cm]  mass_gas(<r)[g]"
        else:
            raise ValueError("unknown profile name: %s" % profile)
        np.savetxt(outfile, data, header=header)


def main():
    parser = argparse.ArgumentParser(
            description="Calculate the mass and potential profiles")
    parser.add_argument("config", nargs="?", default="mass_potential.conf",
                        help="config for mass and potential profiles " +
                             "calculation (default: mass_potential.conf)")
    args = parser.parse_args()

    config = ConfigObj(args.config)

    redshift = config.as_float("redshift")
    pixel = ChandraPixel(z=redshift)

    ne_profile = np.loadtxt(config["ne_profile"])

    cf_profile = np.loadtxt(config["cf_profile"])
    cf_unit = "pixel"
    try:
        cf_unit = config["cf_unit"]
    except KeyError:
        pass
    if cf_unit == "pixel":
        conv_factor = pixel.get_length().to(au.cm).value
    else:
        conv_factor = au.Unit(cf_unit).to(au.cm)
    cf_profile[:, 0] *= conv_factor

    t_profile = np.loadtxt(config["t_profile"])
    t_unit = "pixel"
    try:
        t_unit = config["t_unit"]
    except KeyError:
        pass
    if t_unit == "pixel":
        conv_factor = pixel.get_length().to(au.cm).value
    else:
        conv_factor = au.Unit(t_unit).to(au.cm)
    t_profile[:, 0] *= conv_factor

    density_profile = DensityProfile(density=ne_profile,
                                     density_type="electron")
    density_profile.load_cf_profile(cf_profile)
    density_profile.load_t_profile(t_profile)
    density_profile.calc_gas_density()
    density_profile.interpolate()
    density_profile.gen_radius(num=config.as_int("num_dp"))
    density_profile.calc_mass_gas(verbose=True)
    density_profile.save(profile="mass_gas", outfile=config["m_gas_profile"])


if __name__ == "__main__":
    main()
