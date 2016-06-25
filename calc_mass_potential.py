#!/usr/bin/env python3
#
# Weitian LI
# Created: 2016-06-24
# Updated: 2016-06-25
#
# Change logs:
# 2016-06-25:
#   * Use 'InterpolatedUnivariateSpline' instead of 'interp1d'
#   * Implement 'calc_mass_total()'
#   * Update documentation
# 2016-06-24:
#   * Update method 'gen_radius()'
#

"""
Calculate the (gas and gravitational) mass profile and gravitational
potential profile from the electron number density profile, under the
assumption of hydrostatic equilibrium.

The electron density profile and temperature profile are required.

Assuming that the gas is in hydrostatic equilibrium with the gravitational
potential and a spherically-symmetric distribution of the gas, we can
write the hydrostatic equilibrium equation (HEE) of the ICM as
(ref.[1], eq.(6)):
    derivative(P_gas, r) / rho_gas = - derivative(phi, r)
                                   = - G M_tot(<r) / r^2
where,
    phi: gravitational potential;
    G: gravitational constant;
    rho_gas: gas mass density:
        rho_gas = mu * m_atom * n_gas
    P_gas: gas pressure:
        P_gas = rho_gas * k_B * T_gas / (mu * m_atom) = n_gas * k_B * T_gas
    mu: mean molecular weight in a.m.u (i.e., m_atom) (~ 0.6)
    m_atom: atom mass unit
    n_gas: gas number density; sum of the electron and proton densities
    k_B: Boltzmann constant
    T_gas: gas temperature

Solve the above equation, we can get the total mass of X-ray luminous
galaxy clusters (ref.[1], eq.(7)):
    M_tot(<r) = - (k_B * T_gas(r) * r) / (mu * m_atom * G) *
                (derivative(log(T_gas), log(r)) +
                 derivative(log(n_gas), log(r)))

Note that the second part (the derivatives) is DIMENSIONLESS, since
    d(log(X)) = d(X) / X
Also note that ('R' is a ratio constant):
    d(log(n_gas)) = d(log(R*n_e)) = d(log(n_e))
And and 'log' derivative can be calculated as:
    derivative(log(X(r)), log(r)) = (r / X(r)) * derivative(X(r), r)

Note that 'kT' has dimension of energy.  Therefore, if the gas temperature
is given in 'keV', then the 'kT' should be substitute as a whole.

For example:
    (1.0 keV) * (1.0 kpc) / (0.6 * m_atom * G) ~= 3.7379e10 [ Msun ]
which is consistent with the formula of (ref.[2], eq.(3))


References:
[1] Ettori et al., 2013, Space Science Review, 177, 119-154
[2] Walker et al., 2012, MNRAS, 422, 3503


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
num_dp = 1000

# output gas mass profile
m_gas_profile = mass_gas_profile.txt

# output total (gravitational) mass profile
m_total_profile = mass_total_profile.txt

# output gravitational potential profile
potential_profile = potential_profile.txt
------------------------------------------------------------
"""

import argparse

import numpy as np
import astropy.units as au
import astropy.constants as ac
import scipy.interpolate as interpolate
import scipy.integrate as integrate
from scipy.misc import derivative
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

    def interpolate(self, verbose=False):
        """
        Interpolate the density profile, cooling function profile,
        and temperature profile.

        NOTE:
        * Simple interpolation (`scipy.interpolate.interp1d` of kinds
          `linear` or `quadratic` or `cubic`) may lead to a severely
          oscillating curve (e.g., electron density profile), which
          further cause problems for following mass calculation
          due to the needs to take the derivative.
        * `InterpolatedUnivariateSpline` or `UnivariateSpline` is a
          better choice than the above `interp1d`.
        * Allow cooling function profile and temperature profile to be
          extrapolated by filling with the last value if necessary.
        """
        # density profile
        # insert a data point at radius of zero
        if verbose:
            print("Interpolating density profile ...")
        self.d_interp = interpolate.InterpolatedUnivariateSpline(
            x=np.concatenate([[0.0], self.r]),
            y=np.concatenate([[self.d[0]], self.d]))
        if self.ne is not None:
            if verbose:
                print("Interpolating electron number density profile ...")
            self.ne_interp = interpolate.InterpolatedUnivariateSpline(
                x=np.concatenate([[0.0], self.r]),
                y=np.concatenate([[self.ne[0]], self.ne]))
        if self.rho_gas is not None:
            if verbose:
                print("Interpolating gas mass density profile ...")
            self.rho_gas_interp = interpolate.InterpolatedUnivariateSpline(
                x=np.concatenate([[0.0], self.r]),
                y=np.concatenate([[self.rho_gas[0]], self.rho_gas]))
        # cooling function profile
        if verbose:
            print("Interpolating cooling function profile ...")
        self.cf_interp = interpolate.interp1d(
            x=self.cf_radius, y=self.cf_value, kind="linear",
            bounds_error=False, fill_value=self.cf_value[-1],
            assume_sorted=True)
        # temperature profile
        if verbose:
            print("Interpolating temperature profile ...")
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
            if verbose and (i+1) % 50 == 0:
                print("%d..." % (i+1), end="", flush=True)
            # integrated gas mass [ g ]
            m_gas[i] = integrate.quad(_f_rho_gas, a=0.0, b=r,
                                      epsabs=1e-5*au.kpc.to(au.cm),
                                      epsrel=1e-3)[0]
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
        if self.t_radius is None or self.t_value is None:
            raise ValueError("temperature profile required")
        #
        # calculate the coefficient of the total mass formula
        # M0 = (k_B * T_gas(r) * r) / (mu * m_atom * G)
        M0 = ((1.0*au.keV) * (1.0*au.cm) /
              (AstroParams.mu * ac.u * ac.G)).to(au.g).value
        m_total = np.zeros(self.radius.shape)
        if verbose:
            print("Calculating the total mass profile (#%d): ..." %
                  len(m_total), end="", flush=True)
        for i, r in enumerate(self.radius):
            if verbose and (i+1) % 100 == 0:
                print("%d..." % (i+1), end="", flush=True)
            # enclosed total mass [ g ]
            m_total[i] = - M0 * self.t_interp(r) * r * (
                ((r / self.t_interp(r)) *
                 derivative(self.t_interp, r, dx=0.01*au.kpc.to(au.cm))) +
                ((r / self.ne_interp(r)) *
                 derivative(self.ne_interp, r, dx=0.01*au.kpc.to(au.cm))))
        if verbose:
            print("DONE!", flush=True)
        self.m_total = m_total
        return m_total

    def plot(self, profile, ax=None, fig=None):
        pass

    def save(self, profile, outfile):
        if profile == "mass_gas":
            data = np.column_stack([self.radius,
                                    self.radius_err,
                                    self.m_gas])
            header = "radius[cm]  radius_err[cm]  mass_gas(<r)[g]"
        elif profile == "mass_total":
            data = np.column_stack([self.radius,
                                    self.radius_err,
                                    self.m_total])
            header = "radius[cm]  radius_err[cm]  mass_total(<r)[g]"
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
    density_profile.calc_electron_density()
    density_profile.calc_gas_density()
    density_profile.interpolate(verbose=True)
    density_profile.gen_radius(num=config.as_int("num_dp"))
    density_profile.calc_mass_gas(verbose=True)
    density_profile.save(profile="mass_gas",
                         outfile=config["m_gas_profile"])
    density_profile.calc_mass_total(verbose=True)
    density_profile.save(profile="mass_total",
                         outfile=config["m_total_profile"])


if __name__ == "__main__":
    main()
