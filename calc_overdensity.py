#!/usr/bin/env python3
#
# Aaron LI
# Created: 2016-06-30
# Updated: 2016-06-30
#

"""
Calculate the overdensity profile, and from which to calculate the
R_{500} (defined as the radius of the sphere that encloses a mean total
mass density of 500 times the critical density at the cluster's redshift)
and M_gas_{500}/M_{500} (the enclosed gas/total mass by a sphere of radius
R_{500}).


References:
[1] Ettori et al., 2013, Space Science Review, 177, 119-154


Sample configuration file:
------------------------------------------------------------
## Configuration for `calc_overdensity.py`
## Date: 2016-06-30

# redshift of the source (critical density)
redshift = <REDSHIFT>

# gas mass profile
m_gas_profile = mass_gas_profile.txt

# output total (gravitational) mass profile
m_total_profile = mass_total_profile.txt

# number of times w.r.t the critical density
delta = 500

# output results in JSON format
outfile = overdensity.json

# output overdensity profile
overdensity_profile = overdensity_profile.txt
------------------------------------------------------------
"""

import argparse
import json
from collections import OrderedDict

import numpy as np
import scipy.optimize as optimize
import astropy.units as au
from astropy.cosmology import FlatLambdaCDM
from configobj import ConfigObj

import rpy2.robjects as ro
from rpy2.robjects.packages import importr

from astro_params import AstroParams


class MassProfile:
    """
    Cluster's gas/total integrated mass profile.

    The total/gravitational mass profile is required to calculate
    the overdensity profile, from which the R_{delta} is then determined.
    """
    # supported types of mass profile
    MASS_TYPES = ["total", "gas"]
    # available splines
    SPLINES = ["mass", "overdensity"]
    # input mass data: [r, r_err, m]
    r = None
    r_err = None
    m = None
    # redshift of the object
    redshift = None
    # fitted smoothing spline
    m_spline = None
    m_spline_log10 = None
    od_spline = None
    od_spline_log10 = None
    # call R through `rpy2`
    mgcv = importr("mgcv")

    def __init__(self, mass, mass_type="total"):
        self.load_data(data=mass, mass_type=mass_type)

    def load_data(self, data, mass_type="total"):
        if mass_type not in self.MASS_TYPES:
            raise ValueError("invalid mass_types: %s" % mass_type)
        # 3-column mass profile: [r, r_err, mass]
        self.r = data[:, 0].copy()
        self.r_err = data[:, 1].copy()
        self.m = data[:, 2].copy()
        self.mass_type = mass_type

    def calc_overdensity(self, z, verbose=True):
        """
        Calculate the overdensity profile from the total/gravitational
        mass profile.

        The overdensity is the ratio of the enclosed mean total mass
        density to the critical density at the source's redshift.
        """
        if self.mass_type != "total":
            raise ValueError("total mass profile is required")
        #
        if verbose:
            print("Calculating the overdensity profile ...")
        overdensity = np.zeros(self.r.shape)
        # critical density
        cosmo = FlatLambdaCDM(H0=AstroParams.H0, Om0=AstroParams.OmegaM0)
        d_crit = cosmo.critical_density(z).cgs.value  # [ g/cm^3 ]
        for i, r in enumerate(self.r):
            volume = (4.0/3.0) * np.pi * r**3
            overdensity[i] = self.m[i] / volume / d_crit
        self.overdensity = overdensity
        return overdensity

    def calc_radius_delta(self, delta):
        """
        Calculate the radius at which the overdensity is delta.
        """
        if self.od_spline is None:
            self.fit_spline(spline="overdensity", log10=True)
        if min(self.overdensity) > delta:
            raise ValueError("min(overdensity) > %d" % delta)
        r = optimize.newton(
            lambda x: self.eval_spline("overdensity", x) - delta,
            x0=500.0*au.kpc.to(au.cm),
            tol=1e-2*au.kpc.to(au.cm))
        return r

    def calc_mass_delta(self, r_delta):
        if self.m_spline is None:
            self.fit_spline(spline="mass", log10=True)
        return self.eval_spline(spline="mass", x=r_delta)

    def save(self, outfile):
        """
        Save calculated overdensity profile.
        """
        data = np.column_stack([self.r,
                                self.r_err,
                                self.m])
        header = "radius[cm]  radius_err[cm]  overdensity"
        np.savetxt(outfile, data, header=header)

    def fit_spline(self, spline, log10):
        """
        Fit a smoothing line to the specified spline data,
        by utilizing the R `mgcv::gam()` function.

        If 'log10' is True, the input data are first applied the log-scale
        transform, and then fitted by the smoothing spline.

        The fitted spline allows extrapolation.
        """
        if spline not in self.SPLINES:
            raise ValueError("invalid spline: %s" % spline)
        #
        if spline == "mass":
            # input gas/total mass profile
            if log10:
                x = ro.FloatVector(np.log10(self.r))
                y = ro.FloatVector(np.log10(self.m))
                self.m_spline_log10 = True
            else:
                x = ro.FloatVector(self.r)
                y = ro.FloatVector(self.m)
                self.m_spline_log10 = False
            self.m_spline = self.mgcv.gam(
                ro.Formula("y ~ s(x, bs='ps')"),
                data=ro.DataFrame({"x": x, "y": y}),
                method="REML")
        elif spline == "overdensity":
            # calculated overdensity profile
            if log10:
                x = ro.FloatVector(np.log10(self.r))
                y = ro.FloatVector(np.log10(self.overdensity))
                self.od_spline_log10 = True
            else:
                x = ro.FloatVector(self.radius)
                y = ro.FloatVector(self.rho_total)
                self.od_spline_log10 = False
            self.od_spline = self.mgcv.gam(
                ro.Formula("y ~ s(x, bs='ps')"),
                data=ro.DataFrame({"x": x, "y": y}),
                method="REML")
        else:
            raise ValueError("invalid spline: %s" % spline)

    def eval_spline(self, spline, x):
        """
        Evaluate the specified spline at the supplied positions.
        Also check whether the spline was fitted in the log-scale space,
        and transform the evaluated values if necessary.
        """
        x = np.array(x)
        if x.shape == ():
            x = x.reshape((1,))
        if spline == "mass":
            spl = self.m_spline
            log10 = self.m_spline_log10
        elif spline == "overdensity":
            spl = self.od_spline
            log10 = self.od_spline_log10
        else:
            raise ValueError("invalid spline: %s" % spline)
        #
        if log10:
            x_new = ro.ListVector({"x": ro.FloatVector(np.log10(x))})
            y_pred = self.mgcv.predict_gam(spl, newdata=x_new)
            y_pred = 10 ** np.array(y_pred)
        else:
            x_new = ro.ListVector({"x": ro.FloatVector(x)})
            y_pred = self.mgcv.predict_gam(spl, newdata=x_new)
            y_pred = np.array(y_pred)
        #
        if len(y_pred) == 1:
            return y_pred[0]
        else:
            return y_pred


def main():
    parser = argparse.ArgumentParser(
            description="Calculate overdensity profile and R_{500} etc.")
    parser.add_argument("config", nargs="?", default="overdensity.conf",
                        help="config for overdensity profile and R_{500} " +
                             "etc. calculation (default: overdensity.conf)")
    args = parser.parse_args()

    config = ConfigObj(args.config)
    redshift = config.as_float("redshift")
    m_gas_data = np.loadtxt(config["m_gas_profile"])
    m_total_data = np.loadtxt(config["m_total_profile"])
    delta = list(map(int, config.as_list("delta")))

    m_total_profile = MassProfile(mass=m_total_data, mass_type="total")
    m_total_profile.calc_overdensity(z=redshift, verbose=True)
    m_total_profile.save(outfile=config["overdensity_profile"])

    m_gas_profile = MassProfile(mass=m_gas_data, mass_type="gas")

    results = OrderedDict()
    results["redshift"] = redshift

    for d in delta:
        r_delta = m_total_profile.calc_radius_delta(delta=d)
        m_total_delta = m_total_profile.calc_mass_delta(r_delta)
        m_gas_delta = m_gas_profile.calc_mass_delta(r_delta)
        results["R%d[kpc]" % d] = r_delta * au.cm.to(au.kpc)
        results["Mtotal%d[Msun]" % d] = m_total_delta * au.g.to(au.solMass)
        results["Mgas%d[Msun]" % d] = m_gas_delta * au.g.to(au.solMass)

    results_json = json.dumps(results, indent=2)
    print(results_json)
    open(config["outfile"], "w").write(results_json+"\n")


if __name__ == "__main__":
    main()
