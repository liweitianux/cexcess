#!/usr/bin/env python3
#
# Copyright (c) 2016 Aaron LI
# MIT license
#
# Created: 2016-06-24
#
# Change logs:
# 2016-07-13:
#   * Rename from 'calc_mass_potential.py' to 'calc_mass.py'
#   * Split out the potential profile calculation -> 'calc_potential.py'
# 2016-07-11:
#   * Use a default config to allow a minimal user config
# 2016-07-10:
#   * Allow disable the calculation of potential profile
#   * Use class 'SmoothSpline' from module 'spline.py'
# 2016-07-04:
#   * Remove unnecessary configuration options
#   * Update radii unit to be "kpc", mass unit to be "Msun"
# 2016-06-29:
#   * Update "plot()" to support electron number density profile
# 2016-06-28:
#   * Implement plot function
#   * Adjust integration tolerances and progress report
#   * Fit smoothing splines to profiles using R `mgcv::gam()`
# 2016-06-27:
#   * Implement potential profile calculation:
#     methods 'calc_density_total()' and 'calc_potential()'
# 2016-06-26:
#   * Add document on gravitational potential calculation
# 2016-06-25:
#   * Rename method 'interpolate()' to 'fit_spline()'
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

Note that the second part (i.e., the derivatives) is DIMENSIONLESS, since
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

------------------------------------------------------------
References:
[1] Ettori et al., 2013, Space Science Review, 177, 119-154
[2] Walker et al., 2012, MNRAS, 422, 3503
"""


import argparse

import numpy as np
import astropy.units as au
import astropy.constants as ac
import scipy.integrate as integrate
from scipy.misc import derivative
from configobj import ConfigObj
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from astro_params import AstroParams, ChandraPixel
from projection import Projection
from spline import SmoothSpline

plt.style.use("ggplot")

config_default = """
## Configuration for `calc_mass.py`

# electron density profile
ne_profile = ne_profile.txt

# cooling function profile
cf_profile = coolfunc_profile.txt

# temperature profile
t_profile = t_profile.txt

# number of data points for the output profiles (interpolation)
num_dp = 1000

# output gas mass profile
m_gas_profile = mass_gas_profile.txt
m_gas_profile_image = mass_gas_profile.png

# output total (gravitational) mass profile
m_total_profile = mass_total_profile.txt
m_total_profile_image = mass_total_profile.png

# output total mass density profile
rho_total_profile = rho_total_profile.txt
rho_total_profile_image = rho_total_profile.png

# output gravitational potential profile
# NOTE: to disable potential calculation, do not specified the output files
potential_profile = potential_profile.txt
potential_profile_image = potential_profile.png
"""


class DensityProfile:
    """
    Utilize the 3D (electron number or gas mass) density profile to
    calculate the following quantities:
    * 2D projected surface brightness (requires cooling function profile)
    * gas mass profile (integrated, M_gas(<r))
    * gravitational mass profile (M(<r); requires temperature profile)
    * gravitational potential profile (cut at the largest available radius)

    NOTE:
    * The radii should have unit [ kpc ]
    * The density should have unit [ cm^-3 ] or [ g cm^-3 ]
    """
    # available splines
    SPLINES = ["density", "electron", "rho_gas",
               "cooling_function", "temperature",
               "mass_total", "rho_total"]
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
    # generated radial data points for profile calculation
    radius = None
    radius_err = None
    # gas mass profile: M_gas(<r); same length as the above 'radius'
    m_gas = None
    # total (gravitational) mass profile: M_total(<r)
    m_total = None
    # total mass density profile (required by potential calculation)
    rho_total = None
    # fitted spline to the profiles
    d_spline = None
    ne_spline = None
    rho_gas_spline = None
    cf_spline = None
    t_spline = None
    m_total_spline = None
    rho_total_spline = None

    def __init__(self, density, density_type="electron"):
        self.load_data(data=density, density_type=density_type)

    def load_data(self, data, density_type="electron"):
        if density_type not in self.DENSITY_TYPES:
            raise ValueError("invalid density_types: %s" % density_type)
        # 3-column density profile: r[kpc], r_err[kpc], density
        self.r = data[:, 0].copy()
        self.r_err = data[:, 1].copy()
        self.d = data[:, 2].copy()
        self.density_type = density_type

    def load_cf_profile(self, data):
        if data.shape[1] == 2:
            # 2-column cooling function profile: r[kpc], cf[flux/EM]
            self.cf_radius = data[:, 0].copy()
            self.cf_value = data[:, 1].copy()
        elif data.shape[1] == 3:
            # 3-column cooling function profile: r[kpc], r_err, cf[flux/EM]
            self.cf_radius = data[:, 0].copy()
            self.cf_value = data[:, 2].copy()
        else:
            raise ValueError("invalid cooling function profile data")

    def load_t_profile(self, data):
        if data.shape[1] == 2:
            # 2-column temperature profile: r[kpc], T[keV]
            self.t_radius = data[:, 0].copy()
            self.t_value = data[:, 1].copy()
        elif data.shape[1] == 3:
            # 3-column temperature profile: r[kpc], r_err[kpc], T[keV]
            self.t_radius = data[:, 0].copy()
            self.t_value = data[:, 2].copy()
        else:
            raise ValueError("invalid temperature profile data")

    def calc_density_electron(self):
        """
        Calculate the electron number density from the gas mass density
        if necessary, and fit a smoothing spline to it.
        """
        if self.density_type == "electron":
            self.ne = self.d.copy()
        elif self.density_type == "gas":
            self.ne = self.d / AstroParams.m_atom / AstroParams.mu_e
        # fit a smoothing spline
        self.fit_spline(spline="electron", log10=["x", "y"])
        return self.ne

    def calc_density_gas(self):
        """
        Calculate the gas mass density from the electron number density
        if necessary, and fit a smoothing spline to it.
        """
        if self.density_type == "electron":
            self.rho_gas = self.d * AstroParams.mu_e * AstroParams.m_atom
        elif self.density_type == "gas":
            self.rho_gas = self.d.copy()
        # fit a smoothing spline
        self.fit_spline(spline="rho_gas", log10=["x", "y"])
        return self.rho_gas

    def calc_brightness(self):
        """
        Project the electron number density or gas mass density profile
        to calculate the 2D surface brightness profile.
        """
        if self.cf_radius is None or self.cf_value is None:
            raise ValueError("cooling function profile missing")
        if self.cf_spline is None:
            self.fit_spline(spline="cooling_function", log10=[])
        #
        ne = self.calc_density_electron()
        # flux per unit volume
        cf_new = self.eval_spline(spline="cooling_function", x=self.r)
        flux = cf_new * ne ** 2 / AstroParams.ratio_ne_np
        # project the 3D flux into 2D brightness
        rout = (self.r + self.r_err) * au.kpc.to(au.cm)
        projector = Projection(rout)
        brightness = projector.project(flux)
        return brightness

    def fit_spline(self, spline, log10=[]):
        if spline not in self.SPLINES:
            raise ValueError("invalid spline: %s" % spline)
        #
        if spline == "density":
            # given density profile (either electron / gas mass density)
            x = self.r
            y = self.d
            spl = "d_spline"
        elif spline == "electron":
            # input electron number density profile
            x = self.r
            y = self.ne
            spl = "ne_spline"
        elif spline == "rho_gas":
            # input gas mass density profile
            x = self.r
            y = self.rho_gas
            spl = "rho_gas_spline"
        elif spline == "cooling_function":
            # input cooling function profile
            x = self.cf_radius
            y = self.cf_value
            spl = "cf_spline"
        elif spline == "temperature":
            # input temperature profile
            x = self.t_radius
            y = self.t_value
            spl = "t_spline"
        elif spline == "mass_total":
            # calculated total/gravitational mass profile
            x = self.radius
            y = self.m_total
            spl = "m_total_spline"
        elif spline == "rho_total":
            # calculated total mass density profile
            x = self.radius
            y = self.rho_total
            spl = "rho_total_spline"
        else:
            raise ValueError("invalid spline: %s" % spline)
        setattr(self, spl, SmoothSpline(x=x, y=y))
        getattr(self, spl).fit(log10=log10)

    def eval_spline(self, spline, x):
        """
        Evaluate the specified spline at the supplied positions.
        Also check whether the spline was fitted in the log-scale space,
        and transform the evaluated values if necessary.
        """
        if spline == "density":
            spl = self.d_spline
        elif spline == "electron":
            spl = self.ne_spline
        elif spline == "rho_gas":
            spl = self.rho_gas_spline
        elif spline == "cooling_function":
            spl = self.cf_spline
        elif spline == "temperature":
            spl = self.t_spline
        elif spline == "mass_total":
            spl = self.m_total_spline
        elif spline == "rho_total":
            spl = self.rho_total_spline
        else:
            raise ValueError("invalid spline: %s" % spline)
        return spl.eval(x)

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
            return 4*np.pi * r**2 * self.eval_spline(spline="rho_gas", x=r)
        #
        m_gas = np.zeros(self.radius.shape)
        if verbose:
            print("Calculating the gas mass profile (#%d): ..." %
                  len(self.radius), end="", flush=True)
        c = au.kpc.to(au.cm)**3 * au.g.to(au.solMass)
        for i, r in enumerate(self.radius):
            if verbose and (i+1) % 50 == 0:
                print("%d..." % (i+1), end="", flush=True)
            # enclosed gas mass [ Msun ]
            m_gas[i] = c * integrate.quad(_f_rho_gas, a=0.0, b=r,
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
        if self.t_spline is None:
            self.fit_spline(spline="temperature", log10=[])
        #
        # calculate the coefficient of the total mass formula
        # M0 = (k_B * T_gas(r) * r) / (mu * m_atom * G)
        M0 = ((1.0*au.keV) * (1.0*au.kpc) /
              (AstroParams.mu * ac.u * ac.G)).to(au.solMass).value
        m_total = np.zeros(self.radius.shape)
        if verbose:
            print("Calculating the total mass profile (#%d): ..." %
                  len(self.radius), end="", flush=True)
        for i, r in enumerate(self.radius):
            if verbose and (i+1) % 100 == 0:
                print("%d..." % (i+1), end="", flush=True)
            T = self.eval_spline(spline="temperature", x=r)
            dT_dr = derivative(lambda r: self.eval_spline("temperature", r),
                               r, dx=0.01)
            ne = self.eval_spline(spline="electron", x=r)
            dne_dr = derivative(lambda r: self.eval_spline("electron", r),
                                r, dx=0.01)
            # enclosed total mass [ Msun ]
            m_total[i] = - M0 * T * r * (((r / T) * dT_dr) +
                                         ((r / ne) * dne_dr))
        if verbose:
            print("DONE!", flush=True)
        self.m_total = m_total
        return m_total

    def calc_density_total(self, verbose=True):
        """
        Calculate the total mass density profile, which is required to
        calculate the following gravitational potential profile.
        """
        if self.m_total_spline is None:
            self.fit_spline(spline="mass_total", log10=["x", "y"])
        #
        if verbose:
            print("Calculating the total mass density profile ...")
        rho_total = np.zeros(self.radius.shape)
        # unit conversion: Msun/kpc^3 -> g/cm^3
        c = au.solMass.to(au.g) / au.kpc.to(au.cm)**3
        for i, r in enumerate(self.radius):
            dM_dr = derivative(lambda r: self.eval_spline("mass_total", r),
                               r, dx=0.01)
            rho_total[i] = (dM_dr / (4 * np.pi * r**2)) * c
        self.rho_total = rho_total
        return rho_total

    def plot(self, profile, with_spline=True, ax=None, fig=None):
        x = self.radius
        xlabel = "3D Radius"
        xunit = "kpc"
        xscale = "log"
        yscale = "log"
        x_spl, y_spl = None, None
        if profile == "electron":
            x = self.r
            y = self.ne
            ylabel = "Electron number density"
            yunit = "cm$^{-3}$"
            if with_spline:
                x_spl = self.radius
                y_spl = self.eval_spline(spline="electron", x=self.radius)
        elif profile == "mass_gas":
            y = self.m_gas
            ylabel = "Gas mass"
            yunit = "M$_{\odot}$"
        elif profile == "mass_total":
            y = self.m_total
            ylabel = "Total mass"
            yunit = "M$_{\odot}$"
        elif profile == "rho_total":
            y = self.rho_total
            ylabel = "Total mass density"
            yunit = "g/cm$^3$"
        else:
            raise ValueError("unknown profile name: %s" % profile)
        #
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        ax.plot(x, y, linewidth=2)
        if with_spline and y_spl is not None:
            ax.plot(x_spl, y_spl, linewidth=2, linestyle="dashed")
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        ax.set_xlim(min(x), max(x))
        y_min, y_max = min(y), max(y)
        y_min = y_min/1.2 if y_min > 0 else y_min*1.2
        y_max = y_max*1.2 if y_max > 0 else y_max/1.2
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel(r"%s (%s)" % (xlabel, xunit))
        ax.set_ylabel(r"%s (%s)" % (ylabel, yunit))
        fig.tight_layout()
        return (fig, ax)

    def save(self, profile, outfile):
        if profile == "mass_gas":
            data = np.column_stack([self.radius,
                                    self.radius_err,
                                    self.m_gas])
            header = "radius[kpc]  radius_err[kpc]  mass_gas(<r)[Msun]"
        elif profile == "mass_total":
            data = np.column_stack([self.radius,
                                    self.radius_err,
                                    self.m_total])
            header = "radius[kpc]  radius_err[kpc]  mass_total(<r)[Msun]"
        elif profile == "rho_total":
            data = np.column_stack([self.radius,
                                    self.radius_err,
                                    self.rho_total])
            header = "radius[kpc]  radius_err[kpc]  density_total[g/cm^3]"
        else:
            raise ValueError("unknown profile name: %s" % profile)
        np.savetxt(outfile, data, header=header)


def main():
    parser = argparse.ArgumentParser(
            description="Calculate the gas/total mass profiles")
    parser.add_argument("config", nargs="?",
                        help="additional config")
    args = parser.parse_args()

    config = ConfigObj(config_default.splitlines())
    if args.config is not None:
        config_user = ConfigObj(open(args.config))
        config.merge(config_user)

    ne_profile = np.loadtxt(config["ne_profile"])
    cf_profile = np.loadtxt(config["cf_profile"])
    t_profile = np.loadtxt(config["t_profile"])

    density_profile = DensityProfile(density=ne_profile,
                                     density_type="electron")
    density_profile.load_cf_profile(cf_profile)
    density_profile.load_t_profile(t_profile)
    density_profile.calc_density_electron()
    density_profile.calc_density_gas()
    density_profile.gen_radius(num=config.as_int("num_dp"))

    density_profile.calc_mass_gas(verbose=True)
    density_profile.save(profile="mass_gas",
                         outfile=config["m_gas_profile"])
    fig = Figure(figsize=(10, 8))
    FigureCanvas(fig)
    ax = fig.add_subplot(1, 1, 1)
    density_profile.plot(profile="mass_gas", ax=ax, fig=fig)
    fig.savefig(config["m_gas_profile_image"], dpi=150)

    density_profile.calc_mass_total(verbose=True)
    density_profile.save(profile="mass_total",
                         outfile=config["m_total_profile"])
    fig = Figure(figsize=(10, 8))
    FigureCanvas(fig)
    ax = fig.add_subplot(1, 1, 1)
    density_profile.plot(profile="mass_total", ax=ax, fig=fig)
    fig.savefig(config["m_total_profile_image"], dpi=150)

    density_profile.calc_density_total(verbose=True)
    density_profile.save(profile="rho_total",
                         outfile=config["rho_total_profile"])
    fig = Figure(figsize=(10, 8))
    FigureCanvas(fig)
    ax = fig.add_subplot(1, 1, 1)
    density_profile.plot(profile="rho_total", ax=ax, fig=fig)
    fig.savefig(config["rho_total_profile_image"], dpi=150)


if __name__ == "__main__":
    main()
