#!/usr/bin/env python3
#
# Aaron LI
# Created: 2016-07-04
# Updated: 2016-07-11
#
# Change logs:
# 2016-07-11:
#   * Use a default config to allow a minimal user config
# 2016-07-04:
#   * Set default "rcut=3000" for TemperatureProfile.extrapolate()
#

"""
Fit the deprojected ICM temperature data points with a self-proposed
temperature profile model, i.e., the *wang2012* model:
    T(r) = A * (pow(x,n)+xi*a2) / (pow(x,n)+a2) / pow(1+x*x/a3/a3, beta) + T0

With the fitted temperature profile model, we can interpolate and
extrapolate the temperature profile for later mass profile calculation.
"""


import argparse
import json
from collections import OrderedDict

import numpy as np
import astropy.units as au
import lmfit
from configobj import ConfigObj
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from fitting_models import FittingModel
from astro_params import ChandraPixel

plt.style.use("ggplot")

config_default = """
## Configuration for `fit_tprofile.py`

# redshift of the object (for pixel to distance conversion)
redshift = -1

# input temperature profile data file
t_profile_data = t_profile_data.txt

# cut radius to which stop the extrapolation (unit: kpc)
rcut_extrap = 3000

# number of data points for the output temperature profile
num_dp = 1000

# output json file to save the fitting results
t_profile_json = t_profile.json

# output interpolated and extrapolated temperature profile
t_profile = t_profile.txt
t_profile_image = t_profile.png

[model_params]
  # name = initial, lower, upper, variable (FIXED/False to fix the parameter)
  A    = 5.0,  1.0, 500
  n    = 5.0,  0.1, 10
  xi   = 0.3,  0.1, 1.0
  a2   = 2000, 1.0, 1e+05
  a3   = 1000, 100, 3000
  #beta = 0.5,  0.1, 1.0, FIXED
  beta = 0.5,  0.1, 1.0
  T0   = 2.0,  1.0, 5.0
"""


class Wang2012Model(FittingModel):
    """
    *wang2012* model proposed to fit the ICM temperature profile.
    """
    name = "Wang2012 Temperature Profile Model"
    # model parameters
    params = lmfit.Parameters()
    params.add_many(  # (name, value, vary, min, max, expr)
                    ("A",    5.0,  True, 1.0, 500,   None),
                    ("n",    5.0,  True, 0.1, 10,    None),
                    ("xi",   0.3,  True, 0.1, 1.0,   None),
                    ("a2",   2000, True, 1.0, 1.0e5, None),
                    ("a3",   1000, True, 100, 3000,  None),
                    ("beta", 0.5,  True, 0.1, 1.0,   None),
                    ("T0",   2.0,  True, 1.0, 5.0,   None))

    def __init__(self, fit_method="lbfgsb", params=None):
        super().__init__(fit_method=fit_method, params=params, scale=False)

    @staticmethod
    def model(x, params):
        parvals = params.valuesdict()
        A = parvals["A"]
        n = parvals["n"]
        xi = parvals["xi"]
        a2 = parvals["a2"]
        a3 = parvals["a3"]
        beta = parvals["beta"]
        T0 = parvals["T0"]
        return (A * (x**n + xi*a2) / (x**n + a2) /
                ((1 + (x/a3)**2) ** beta) + T0)


class TemperatureProfile:
    """
    Fit the deprojected ICM temperature data points with a temperature
    profile model, and output the interpolated and extrapolated temperature
    profile for later mass profile calculation.

    The input radii have unit "pixel", which are first converted to
    "kpc" and then fitted with the model.

    The output temperature profile also has unit "kpc" for radii.
    """
    # input temperature profile data: [r, r_err, t, t_err]
    r = None
    r_err = None
    t = None
    t_err = None
    # redshift of the source
    z = None
    # `ChandraPixel` instance for unit conversion
    pixel = None
    # flag to indicate whether the units are converted
    units_converted = False
    # model to be fitted
    model = None

    def __init__(self, data, z):
        self.load_data(data)
        self.z = z
        self.pixel = ChandraPixel(z)
        self.model = Wang2012Model()

    def load_data(self, data):
        # 4-column t profile: [r, r_err, temperature, temperature_err]
        self.r = data[:, 0].copy()
        self.r_err = data[:, 1].copy()
        self.t = data[:, 2].copy()
        self.t_err = data[:, 3].copy()

    def convert_units(self):
        """
        Convert the units of input data:
           radius: pixel -> kpc
        """
        if not self.units_converted:
            kpc_per_pixel = self.pixel.get_length().to(au.kpc).value
            self.r *= kpc_per_pixel
            self.r_err *= kpc_per_pixel
            self.units_converted = True

    def fit(self):
        self.model.load_data(xdata=self.r, ydata=self.t,
                             xerr=self.r_err, yerr=self.t_err)
        self.model.fit()

    def extrapolate(self, rcut=3000, num=1000):
        """
        Interpolate and extrapolate the fitted temperature profile.

        The output radii are generated to be linear-evenly distributed.
        """
        self.rcut_extrap = rcut
        self.num_dp = num
        radius = np.linspace(0.0, rcut, num+1)
        rin = radius[:-1]
        rout = radius[1:]
        self.r_extrapolated = (rout + rin) / 2.0
        self.r_err_extrapolated = (rout - rin) / 2.0
        self.t_extrapolated = self.model.f(self.r_extrapolated)

    def report(self, outfile=None):
        """
        Report the temperature profile model fitting results.
        """
        results = OrderedDict([
            ("redshift",    self.z),
            ("rcut_extrap", self.rcut_extrap),
            ("num_dp",      self.num_dp),
            ("model",       self.model.name),
            ("fitting",     self.model.report(rtype="fitting")),
            ("params",      self.model.report(rtype="parameters")),
        ])
        results_json = json.dumps(results, indent=2)
        if outfile is None:
            print(results_json)
        else:
            open(outfile, "w").write(results_json+"\n")

    def save(self, outfile):
        data = np.column_stack([self.r_extrapolated,
                                self.r_err_extrapolated,
                                self.t_extrapolated])
        header = "radius[kpc]  radius_err[kpc]  temperature[keV]"
        np.savetxt(outfile, data, header=header)

    def plot(self, ax=None, fig=None):
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        ax.errorbar(self.r, self.t,
                    xerr=self.r_err, yerr=self.t_err,
                    fmt="none", elinewidth=2, capthick=2)
        # fitted model
        ax.plot(self.r_extrapolated, self.t_extrapolated,
                color="black", linestyle="solid", linewidth=2)
        ax.set_xlabel("Radius (kpc)")
        ax.set_ylabel("Temperature (keV)")
        fig.tight_layout()
        return (fig, ax)


def main():
    parser = argparse.ArgumentParser(
        description="temperature profile fit, interpolate and extrapolate")
    parser.add_argument("config", nargs="?", default="tprofile.conf",
                        help="configuration (default: tprofile.conf")
    args = parser.parse_args()

    config = ConfigObj(config_default.splitlines())
    config_user = ConfigObj(args.config)
    config.merge(config_user)

    tprofile_data = np.loadtxt(config["t_profile_data"])
    redshift = config.as_float("redshift")

    tprofile = TemperatureProfile(tprofile_data, redshift)
    tprofile.convert_units()
    # Load parameters settings from config
    params = config["model_params"]
    for p, value in params.items():
        variable = True
        if len(value) == 4 and value[3].upper() in ["FIXED", "FALSE"]:
            variable = False
        tprofile.model.set_param(name=p, value=float(value[0]),
                                 min=float(value[1]), max=float(value[2]),
                                 vary=variable)
    tprofile.fit()
    tprofile.extrapolate(rcut=config.as_float("rcut_extrap"),
                         num=config.as_int("num_dp"))
    tprofile.report(outfile=config["t_profile_json"])
    tprofile.save(outfile=config["t_profile"])

    fig = Figure(figsize=(10, 8))
    FigureCanvas(fig)
    ax = fig.add_subplot(1, 1, 1)
    tprofile.plot(ax=ax, fig=fig)
    fig.savefig(config["t_profile_image"], dpi=150)


if __name__ == "__main__":
    main()
