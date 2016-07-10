#!/usr/bin/env python3
#
# Aaron LI
# Created: 2016-06-10
# Updated: 2016-07-10
#
# Change logs:
# 2016-07-10:
#   * Use class 'SmoothSpline' from module
# 2016-07-04:
#   * Use model's "report()" method
#   * Add config "sbpexp_rcut"
#   * Rename config "sbpexp_rcut*" to "sbpexp_rignore*"
#   * Save profile radii in unit "kpc"
#   * Update to that cooling function profile's radius in unit "kpc"
# 2016-06-27:
#   * Minor cleanups
#   * Remove obsolete class "DeprojectSBP"
#   * Fit smoothing spline to SBP and cooling function profiles by
#     calling the R `mgcv::gam()`: "fit_spline()" and "eval_spline()"
#   * Update "plot()" to also plot the fitted smoothing spline
# 2016-06-26:
#   * Split out method "save()" for class "BrightnessProfile"
#   * Split classes 'FitModel', 'ABModel' and 'PLCModel' into separate
#     module 'fitting_models.py'
# 2016-06-25:
#   * Use 'InterpolatedUnivariateSpline' instead of 'interp1d'
# 2016-06-24:
#   * Move class 'ChandraPixel' to module 'astro_params.py'
#   * Split class 'Projection' to a separate module 'projection.py'
#   * Move class 'DensityProfile' to tool 'calc_mass_potential.py'
#   * Split class 'AstroParams' to separate module 'astro_params.py'
# 2016-06-23:
#   * Add configuration parameter 'sbpexp_rcut'
#   * Allow extrapolate the cooling function profile
#   * Add plot function to class 'BrightnessProfile'
#   * Update sample configuration file
#   * Remove obsolete class 'SurfaceBrightnessProfile'
# 2016-06-22:
#   * Add class 'DensityProfile', the inversion to 'BrightnessProfile'
#   * Add classes 'AstroParams' and 'BrightnessProfile'
#   * Add class 'ChandraPixel'
#   * Update documentation
# 2016-06-21:
#   * Add document about the gas density derivation
# 2016-06-20:
#   * Use configuration file instead of the tedious command line arguments
# 2016-06-16:
#   * Add methods 'save()', 'report()' and 'plot()' to class 'SBP'
# 2016-06-15:
#   * Add command line arguments
#   * Add class 'SBP' for SBP background subtraction and extrapolation
# 2016-06-14:
#   * Add class 'PLCModel' based on 'FitModel'
#   * Split class 'FitModel' from 'ABModel'
# 2016-06-13:
#   * Add class 'ABModel' to support data scaling
#   * Implement primitive SBP deprojection approach for class 'DeprojectSBP'
#

"""
Deproject the 2D surface brightness profile (SBP) into ???

The SBP deprojection is performed using a non-parametric approach with
regularization which add the constraint that the 3D gas density profile
should be smooth.

=======================================================================

Surface brightness (`SUR_FLUX` column of the dmextract'ed radial profile):
    Brightness: [ photon s^-1 cm^-2 pixel^-2 ]
where the 'cm^-2' is due to the instrumental effective area, and the
'pixel^-2' is corresponding to the solid angle with respect to the source
(i.e., [ arcsec^-2 ]).

The flux has dimension:
    Flux: [ photon s^-1 cm^-2 ]
therefore, the dimension of brightness can also be expressed as:
    Brightness: [ Flux pixel^-2 ] = [ Flux sr^-1 ]

The instrument and (time-normalized) exposure map has dimension:
    [ count photon^-1 cm^2 ]
which is used to convert the instrument-specific counts image into physical-
meaningful flux unit.

Emission measure:
    EM = \int n_e n_H dV ~= (n_e^2 / ratio_eH) V  [ cm^-3 ]
where 'ratio_eH' is the ratio of electron density to proton density (n_H).

APEC normalization returned by XSPEC is simply the *emission measure* of
the gas scaled by the distance:
    eta = (\int n_e n_H dV) / (4 pi (D_A (1+z))^2)
assuming (ref. [4]):
    n_H ~= 0.826 n_e
then the gas density (n_H or n_e) can be calculated.

The flux calculated with the XSPEC `flux` command has dimension:
    Flux: [ photon s^-1 cm^-2 ] or [ erg s^-2 cm^-2 ]

When use XSPEC APEC model to calculate the cooling function (Lambda),
its normalization is calculated with EM = 1, therefore:
    norm = 1e-14 / (4 pi (D_A (1+z))^2) * EM
         = 1e-14 / (4 pi (D_A (1+z))^2)  [ cm^-5 ]
where the 'D_A' is the angular diameter distance which can be simply
calculated from its redshift.

With the Galactic absorption (nH), temperature (varies with radius), and
abundance (assumed constant) been set, the cooling function is derived by
using the XSPEC `flux` command.
Therefore, cooling function has dimension:
    Lambda: [ Flux EM^-1 ]

By deprojecting the surface brightness, the flux per volume can be derived,
and EM can be further obtained by incorporating the cooling function,
and finally the (3D) gas density can be determined.
    (projection): EM * Lambda / A  ->  Brightness
where 'A' is the solid angle (i.e., area covered by the source).
=======================================================================

References:
[1] Croston et al. 2006, A&A, 459, 1007-1019
[2] McLaughlin, 1999, ApJ, 117, 2398-2427
[3] Bouchet, 1995, A&AS, 113, 167
[4] Ettori et al, 2013, Space Science Review, 177, 119-154
[5] AtomDB / APEC model:
    * http://www.atomdb.org/faq.php#DensityXSPECnorm
    * https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/XSmodelApec.html


Sample configuration file:
------------------------------------------------------------
## Configuration for `deproject_sbp.py`
## Date: 2016-06-23

# config file for SBP fitting (e.g., sbpfit.conf)
sbpfit_config = sbpfit.conf

# input cooling function profile
coolfunc_profile = coolfunc_profile.txt

# redshift of the object (for pixel to distance conversion)
redshift = <REDSHIFT>

## SBP extrapolation
# ignorance radius from which the SBP is fitted for extrapolation,
# specified by the ratio w.r.t sbpfit rc (default: 1.2 * rc)
sbpexp_rignore_ratio = 1.2
# or directly specify the ignorance radius (override above) (unit: pixel)
sbpexp_rignore = <RIGNORE>
# cut radius to which stop the extrapolation (unit: kpc)
sbpexp_rcut = <RCUT>
# output of the extrapolated SBP
sbpexp_outfile = sbpexp.csv
# extrapolation model information
sbpexp_json = sbpexp.json
# plot of the SBP extrapolation
sbpexp_image = sbpexp.png

## Density profiles
# deprojected 3D electron number density profile
ne_profile = ne_profile.txt
# deprojected 3D gas mass density profile
rho_gas_profile = rho_gas_profile.txt
# image of the density profiles (electron density and/or gas density)
density_profile_image = density_profile.png
------------------------------------------------------------
"""

import argparse
import json
from collections import OrderedDict

import astropy.units as au
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from configobj import ConfigObj

from astro_params import AstroParams, ChandraPixel
from projection import Projection
from fitting_models import PLCModel
from spline import SmoothSpline

plt.style.use("ggplot")


class SBP:
    """
    X-ray surface brightness profile class.

    This class deals with SBP background subtraction and SBP extrapolation.
    """
    # input SBP data: [r, r_err, s, s_err]
    r = None
    r_err = None
    s = None
    s_err = None
    # uniform background been subtracted
    bkg = None
    # ignorance/minimal radius from which the SBP is fitted to the PLCModel
    rignore = None
    # cut radius where the extrapolation stops
    rcut = None
    # PLCModel instance used to extrapolate the SBP
    plcmodel = None

    def __init__(self, r, r_err=None, s=None, s_err=None, rignore=None):
        self.load_data(r=r, r_err=r_err, s=s, s_err=s_err, rignore=rignore)
        self.plcmodel = PLCModel(scale=True)

    def load_data(self, r, r_err=None, s=None, s_err=None, rignore=None):
        if r.ndim == 2 and r.shape[1] == 4:
            # 4-column data
            self.r = r[:, 0].copy()
            self.r_err = r[:, 1].copy()
            self.s = r[:, 2].copy()
            self.s_err = r[:, 3].copy()
        else:
            self.r = np.array(r)
            self.r_err = np.array(r_err)
            self.s = np.array(s)
            self.s_err = np.array(s_err)
        self.rignore = rignore

    def subtract_bkg(self, bkg):
        """
        Subtract the uniform background from the brightness.

        The value of background can be acquired by fitting the whole
        or core-exclude SBP with model consisting of a plain beta model
        and a constant.
        The "AB model" maybe also applicable.
        """
        self.bkg = bkg
        self.s -= bkg
        self.bkg_subtracted = True

    def extrapolate(self, rignore=None, rcut=None):
        """
        Extrapolate the SBP by assuming that the outer SBP follows the
        following relation:
            S_X = A * r^{-alpha},
        which can be determined by model fitting.

        The SBP is extrapolated to the region where the brightness is
        lower than the current observed minimal brightness by one order
        of magnitude, and the extrapolated SBP bins have the same width
        and relative errors as the last SBP bin observed.

        If the 'rcut' is specified, then the SBP extrapolation stops
        when exceeds that radius.

        Note that the uniform background should be subtracted first.

        Return:
          * self.r_extrapolated
          * self.r_err_extrapolated
          * self.s_extrapolated
          * self.s_err_extrapolated
          * self.mask_extrapolated
        """
        if rignore is not None:
            self.rignore = rignore
        if rcut is not None:
            self.rcut = rcut
        self.mask = self.r >= self.rignore
        self.plcmodel.load_data(xdata=self.r[self.mask],
                                ydata=self.s[self.mask],
                                xerr=self.r_err[self.mask],
                                yerr=self.s_err[self.mask],
                                update_params=True)
        self.plcmodel.set_param("bkg", value=0.0, vary=False)
        self.plcmodel.fit()
        last_r_err = self.r_err[-1]
        last_s = self.s[-1]
        last_s_err = self.s_err[-1]
        #
        r_exp = self.r.tolist()
        r_err_exp = self.r_err.tolist()
        s_exp = self.s.tolist()
        s_err_exp = self.s_err.tolist()
        mask_exp = [False] * len(r_exp)
        # do extrapolation
        r_tmp = r_exp[-1] + 2*r_err_exp[-1]
        s_tmp = self.plcmodel.f(r_tmp)
        while True:
            if rcut is not None and r_tmp > rcut:
                break
            if rcut is None and (s_tmp < last_s / 10.0):
                break
            r_exp.append(r_tmp)
            r_err_exp.append(last_r_err)
            s_exp.append(s_tmp)
            s_err_exp.append(s_tmp * last_s_err / last_s)
            mask_exp.append(True)
            r_tmp = r_exp[-1] + 2*r_err_exp[-1]
            s_tmp = self.plcmodel.f(r_tmp)
        # convert to numpy array
        self.r_extrapolated = np.array(r_exp)
        self.r_err_extrapolated = np.array(r_err_exp)
        self.s_extrapolated = np.array(s_exp)
        self.s_err_extrapolated = np.array(s_err_exp)
        self.mask_extrapolated = np.array(mask_exp)

    def report(self, outfile=None):
        """
        Report the extrapolation model fitting results.
        """
        results = OrderedDict([
            ("bkg",            self.bkg),
            ("bkg_subtracted", self.bkg_subtracted),
            ("rignore",        self.rignore),
            ("rcut",           self.rcut),
            ("model",          self.plcmodel.name),
            ("fitting",        self.plcmodel.report(rtype="fitting")),
            ("params",         self.plcmodel.report(rtype="parameters")),
        ])
        results_json = json.dumps(results, indent=2)
        if outfile is None:
            print(results_json)
        else:
            open(outfile, "w").write(results_json+"\n")

    def plot(self, ax=None, fig=None):
        """
        Make a plot of the SBP extrapolation.
        """
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        # ignored data points
        mask_ignore = np.logical_not(self.mask)
        if np.sum(mask_ignore) > 0:
            ax.errorbar(self.r[mask_ignore], self.s[mask_ignore],
                        xerr=self.r_err[mask_ignore],
                        yerr=self.s_err[mask_ignore],
                        fmt="none", elinewidth=1, capthick=1)
        # data points used to fit the PLC model
        ax.errorbar(self.r[self.mask], self.s[self.mask],
                    xerr=self.r_err[self.mask],
                    yerr=self.s_err[self.mask],
                    fmt="none", elinewidth=2, capthick=2)
        # extrapolated data points
        ax.errorbar(self.r_extrapolated[self.mask_extrapolated],
                    self.s_extrapolated[self.mask_extrapolated],
                    xerr=self.r_err_extrapolated[self.mask_extrapolated],
                    yerr=self.s_err_extrapolated[self.mask_extrapolated],
                    fmt="none", elinewidth=1, capthick=1)
        # original data points without background subtraction
        eb = ax.errorbar(self.r, self.s+self.bkg,
                         xerr=self.r_err, yerr=self.s_err,
                         fmt="none", elinewidth=1, capthick=1)
        eb[-1][0].set_linestyle("dashdot")
        eb[-1][1].set_linestyle("dashdot")
        # PLC model
        mask_fit = self.mask_extrapolated.copy()
        mask_fit[:len(self.mask)] = self.mask
        r_fit = self.r_extrapolated[mask_fit]
        s_fit = self.plcmodel.f(r_fit)
        ax.plot(r_fit, s_fit, color="black", linestyle="solid")
        # adjust layout
        r_min = 1.0
        r_max = self.r_extrapolated[-1] + self.r_err_extrapolated[-1]
        s_min = min(self.s_extrapolated) / 1.2
        s_max = max(self.s_extrapolated + self.s_err_extrapolated) * 1.2
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(r_min, r_max)
        ax.set_ylim(s_min, s_max)
        # labels
        ax.set_xlabel("Radius (%s)" % "pixel")
        ax.set_ylabel(r"Surface Brightness (photons/cm$^2$/pixel$^2$/s)")
        ax.text(x=r_max/1.2, y=s_max/1.2,
                s=r"reduced $\chi^2$: %.2f / %.2f = %.2f" % (
                    self.plcmodel.fitted.chisqr, self.plcmodel.fitted.nfree,
                    self.plcmodel.fitted.chisqr/self.plcmodel.fitted.nfree),
                horizontalalignment="right", verticalalignment="top")
        fig.tight_layout()
        return (fig, ax)

    def get_data(self):
        """
        Get the extrapolated data, for following use.
        """
        return np.column_stack([self.r_extrapolated,
                                self.r_err_extrapolated,
                                self.s_extrapolated,
                                self.s_err_extrapolated])

    def save(self, outfile):
        """
        Save the (extrapolated) SBP to the given output file in CSV format.
        """
        df = pd.DataFrame()
        df["radius"] = self.r_extrapolated
        df["radius_err"] = self.r_err_extrapolated
        df["brightness"] = self.s_extrapolated
        df["brightness_err"] = self.s_err_extrapolated
        df["flag_extrapolation"] = self.mask_extrapolated
        flag_fit = np.zeros(self.mask_extrapolated.shape, dtype=bool)
        flag_fit[:len(self.mask)] = self.mask
        df["flag_fit"] = flag_fit
        df.to_csv(outfile, header=True, index=False)


class BrightnessProfile:
    """
    Calculate the electron number density and/or gas mass density profile
    by deprojecting the observed X-ray surface brightness profile and
    incorporating the cooling function profile.

    NOTE:
    * The radii should have unit [ pixel ] (Chandra pixel)
    * The brightness should have unit [ photon s^-1 cm^-2 pixel^-2 ],
      i.e., [ Flux pixel^-2 ] (radial profile column `SUR_FLUX`)
    """
    # available splines
    SPLINES = ["brightness", "cooling_function"]
    # allowed density profile types
    DENSITY_TYPES = ["electron", "gas"]
    # input SBP data: [r, r_err, s, s_err]
    r = None
    r_err = None
    s = None
    s_err = None
    # redshift of the source
    z = None
    # `ChandraPixel` instance for unit conversion
    pixel = None
    # flag to indicate whether the units are converted
    units_converted = False
    # calculated electron density profile
    ne = None
    # calculated gas mass density profile
    rho_gas = None
    # fitted smoothing spline to the SBP
    s_spline = None
    # fitted smoothing spline to the cooling function profile
    cf_spline = None

    def __init__(self, sbp_data, cf_data, z):
        self.load_data(data=sbp_data)
        self.load_cf_data(data=cf_data)
        self.z = z
        self.pixel = ChandraPixel(z)

    def load_data(self, data):
        # 4-column SBP: [r, r_err, brightness, brightness_err]
        self.r = data[:, 0].copy()
        self.r_err = data[:, 1].copy()
        self.s = data[:, 2].copy()
        self.s_err = data[:, 3].copy()

    def load_cf_data(self, data):
        # 2-column cooling function profile
        self.cf_radius = data[:, 0].copy()
        self.cf_value = data[:, 1].copy()

    def convert_units(self):
        """
        Convert the units of SBP:
           radius: pixel -> cm
           brightness: Flux / pixel**2 -> Flux / cm**2

        Convert the units of cooling function profile:
            radius: kpc -> cm
        """
        if not self.units_converted:
            cm_per_pixel = self.pixel.get_length().to(au.cm).value
            self.r *= cm_per_pixel
            self.r_err *= cm_per_pixel
            self.s /= cm_per_pixel**2
            self.s_err /= cm_per_pixel**2
            # cooling function profile: kpc -> cm
            self.cf_radius *= au.kpc.to(au.cm)
            self.units_converted = True

    def get_radius(self):
        return (self.r.copy(), self.r_err.copy())

    def fit_spline(self, spline, log10=[]):
        if spline not in self.SPLINES:
            raise ValueError("invalid spline: %s" % spline)
        #
        if spline == "brightness":
            x = self.r
            y = self.s
            weights = self.s / self.s_err
            spl = "s_spline"
        elif spline == "cooling_function":
            x = self.cf_radius
            y = self.cf_value
            weights = None
            spl = "cf_spline"
        setattr(self, spl, SmoothSpline(x=x, y=y, weights=weights))
        getattr(self, spl).fit(log10=log10)

    def eval_spline(self, spline, x):
        """
        Evaluate the specified spline at the supplied positions.
        Also check whether the spline was fitted in the log-scale space,
        and transform the evaluated values if necessary.
        """
        if spline == "brightness":
            spl = self.s_spline
        elif spline == "cooling_function":
            spl = self.cf_spline
        else:
            raise ValueError("invalid spline: %s" % spline)
        return spl.eval(x)

    def calc_electron_density(self):
        """
        Deproject the surface brightness profile to derive the 3D
        electron number density (and then gas mass density) profile
        by incorporating the cooling function profile.

        unit: [ cm^-3 ] if the units converted for input data
        """
        if self.s_spline is None:
            self.fit_spline(spline="brightness", log10=["x", "y"])
        if self.cf_spline is None:
            self.fit_spline(spline="cooling_function", log10=[])
        #
        s_new = self.eval_spline(spline="brightness", x=self.r)
        cf_new = self.eval_spline(spline="cooling_function", x=self.r)
        #
        projector = Projection(rout=self.r+self.r_err)
        s_deproj = projector.deproject(s_new)
        # emission measure per unit volume
        em_v = s_deproj / cf_new
        ne = np.sqrt(em_v * AstroParams.ratio_ne_np)
        self.ne = ne
        return ne

    def calc_gas_density(self):
        """
        Calculate the gas mass density based the calculated electron
        number density.

        unit: [ g cm^-3 ] if the units converted for input data
        """
        ne = self.calc_electron_density()
        rho = ne * AstroParams.mu_e * AstroParams.m_atom
        self.rho_gas = rho
        return rho

    def save(self, density_type, outfile):
        if density_type == "electron":
            data = np.column_stack([self.r * au.cm.to(au.kpc),
                                    self.r_err * au.cm.to(au.kpc),
                                    self.ne])
            header = "radius[kpc]  radius_err[kpc]  " + \
                     "electron_number_density[cm^-3]"
        elif density_type == "gas":
            data = np.column_stack([self.r * au.cm.to(au.kpc),
                                    self.r_err * au.cm.to(au.kpc),
                                    self.rho_gas])
            header = "radius[kpc]  radius_err[kpc]  " + \
                     "gas_mass_density[g/cm^3]"
        else:
            raise ValueError("unknown density_type: %s" % density_type)
        np.savetxt(outfile, data, header=header)

    def plot(self, ax=None, fig=None, density_type="electron"):
        if density_type not in self.DENSITY_TYPES:
            raise ValueError("invalid density_types: %s" % density_type)
        if density_type == "electron":
            density = self.ne
            d_name = "Deprojected electron number density"
            d_unit = "cm$^{-3}$"
        else:
            density = self.rho_gas
            d_name = "Deprojected gas mass density"
            d_unit = "g/cm$^3$"
        #
        if self.units_converted:
            # convert from [cm] to [kpc]
            r = self.r * au.cm.to(au.kpc)
            r_err = self.r_err * au.cm.to(au.kpc)
            r_unit = "kpc"
            s_unit = "flux/cm$^2$"
        else:
            r = self.r
            r_err = self.r_err
            r_unit = "pixel"
            s_unit = "flux/pixel$^2$"
        #
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        # SBP data points
        eb = ax.errorbar(r, self.s, xerr=r_err, yerr=self.s_err,
                         fmt="none", elinewidth=2, capthick=2,
                         label="Brightness profile")
        # fitted smoothing spline to SBP
        s_new = self.eval_spline(spline="brightness", x=self.r)
        line1, = ax.plot(r, s_new, linestyle="dashed", linewidth=2,
                         label="SBP smoothing spline")
        #
        r_min = 1.0
        r_max = max(r + r_err)
        s_min = min(self.s) / 1.2
        s_max = max(self.s + self.s_err) * 1.2
        ax.set_xlim(r_min, r_max)
        ax.set_ylim(s_min, s_max)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Radius (%s)" % r_unit)
        ax.set_ylabel(r"Surface brightness (%s)" % s_unit)
        # deprojected density profile
        ax2 = ax.twinx()
        line2, = ax2.plot(r, density, color="black",
                          linestyle="solid", linewidth=2,
                          label="Density profile")
        d_min = min(density) / 1.2
        d_max = max(density) * 1.2
        ax2.set_xlim(r_min, r_max)
        ax2.set_ylim(d_min, d_max)
        ax2.set_yscale(ax.get_yscale())
        ax2.set_ylabel(r"%s (%s)" % (d_name, d_unit))
        # legend
        handles = [eb, line1, line2]
        labels = [h.get_label() for h in handles]
        ax.legend(handles=handles, labels=labels, loc=0)
        fig.tight_layout()
        return (fig, ax, ax2)

######################################################################


def main():
    parser = argparse.ArgumentParser(
            description="Deproject the surface brightness profile (SBP)")
    parser.add_argument("config", nargs="?", default="sbpdeproj.conf",
                        help="config for SBP deprojection " +
                             "(default: sbpdeproj.conf)")
    args = parser.parse_args()

    config = ConfigObj(args.config)
    sbpfit_conf = ConfigObj(config["sbpfit_config"])
    try:
        sbpfit_outfile = sbpfit_conf[sbpfit_conf["model"]]["outfile"]
    except KeyError:
        sbpfit_outfile = sbpfit_conf["outfile"]
    sbpfit_results = json.load(open(sbpfit_outfile))
    sbpdata = np.loadtxt(sbpfit_conf["sbpfile"])
    rc = sbpfit_results["params"]["rc"][0]
    bkg = sbpfit_results["params"]["bkg"][0]

    redshift = config.as_float("redshift")
    pixel = ChandraPixel(redshift)

    print("SBP background subtraction and extrapolation ...")
    sbp = SBP(sbpdata)
    # ignorance radius
    rignore = rc * config.as_float("sbpexp_rignore_ratio")
    try:
        rignore = config.as_float("sbpexp_rignore")
    except KeyError:
        pass
    # cut radius where extrapolation stops (unit: kpc)
    try:
        rcut = config.as_float("sbpexp_rcut")
        # convert unit "kpc" -> "pixel"
        rcut /= pixel.get_length().to(au.kpc).value
    except KeyError:
        rcut = None
    sbp.subtract_bkg(bkg=bkg)
    sbp.extrapolate(rignore=rignore, rcut=rcut)
    sbp.save(outfile=config["sbpexp_outfile"])
    sbp.report(outfile=config["sbpexp_json"])
    fig = Figure(figsize=(10, 8))
    FigureCanvas(fig)
    ax = fig.add_subplot(1, 1, 1)
    sbp.plot(ax=ax, fig=fig)
    fig.savefig(config["sbpexp_image"], dpi=150)

    print("SBP deprojection -> density profile ...")
    cf_data = np.loadtxt(config["coolfunc_profile"])
    sbpdata_extrapolated = sbp.get_data()
    brightness_profile = BrightnessProfile(sbp_data=sbpdata_extrapolated,
                                           cf_data=cf_data,
                                           z=redshift)
    brightness_profile.convert_units()
    brightness_profile.calc_electron_density()
    brightness_profile.save(density_type="electron",
                            outfile=config["ne_profile"])
    brightness_profile.calc_gas_density()
    brightness_profile.save(density_type="gas",
                            outfile=config["rho_gas_profile"])
    fig = Figure(figsize=(10, 8))
    FigureCanvas(fig)
    ax = fig.add_subplot(1, 1, 1)
    brightness_profile.plot(ax=ax, fig=fig)
    fig.savefig(config["density_profile_image"], dpi=150)


if __name__ == "__main__":
    main()
