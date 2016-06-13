#!/usr/bin/env python3
#
# Deproject the 2D surface brightness profile (SBP) into the 3D emission
# measure (EM) profile, using the non-parametric approach.
# And the 3D gas density profile can be further derived through the 3D
# EM profile by taking into account the variation of the cooling function
# Lambda(T, Z) with radius.
#
#
# References:
# [1] Croston et al. 2006, A&A, 459, 1007-1019
# [2] McLaughlin, 1999, ApJ, 117, 2398-2427
# [3] Bouchet, 1995, A&AS, 113, 167
#
#
# Weitian LI
# Created: 2016-06-10
# Updated: 2016-06-13
#
# Change logs:
# 2016-06-13:
#   * Add class 'ABModel' to support data scaling
#   * Implement primitive SBP deprojection approach for class 'DeprojectSBP'
#


import argparse

import numpy as np
import scipy.optimize
import lmfit


class Projection:
    """
    Class that deals with projection from 3D volume density to 2D
    surface density and vice versa.

    The inner-most shell/cylinder is assumed to at the center with inner
    radius of ZERO.
    """
    # number of shells/cylinders
    N = 0
    # inner and outer radii of each spherical shell or cylinders
    rin = None
    rout = None
    # projection matrix from 3D volume density to 2D surface density
    proj_mat = None

    def __init__(self, rout):
        self.N = len(rout)
        self.rout = np.array(rout, dtype=float)
        self.rin = np.concatenate([[0.0], self.rout[:-1]])
        self.calc_projection_matrix()

    def __str__(self):
        return "%s: #%d shells: Rout(%s)" % (self.__class__.__name__,
                                             self.N, self.rout)

    def calc_projection_matrix(self):
        """
        Calculate the projection matrix according to the given outer radii.

        Arguments:
        * rout: (vector) outer radius of each SB annulus or spherical shell

        Return:
        * proj_mat: (matrix) an upper triangular matrix with element
                    [i, j] indicate the fraction of the emission from
                    shell j that is observed in annulus i.

        N(R_{i-1}, R_i) * \pi * (R^2_i - R^2_{i-1}) =
            \sum_{j=i}^{m} (n(R_{j-1}, R_j) *
                            V_int(R_{j-1}, R_j; R_{i-1}, R_i))

        References:
        * ref.[1], eq.(1)
        * ref.[2], eq.(A2)
        """
        proj_mat = np.zeros((self.N, self.N))
        for i in range(self.N):
            # loop over each annulus
            rin = self.rin[i]
            rout = self.rout[i]
            area = np.pi * (rout**2 - rin**2)
            for j in range(i, self.N):
                # calculate the contribution from each shell to annulus i
                rin2 = self.rin[j]
                rout2 = self.rout[j]
                v_int = self.intersection_volume(rin2, rout2, rin, rout)
                proj_mat[i, j] = v_int / area
        self.proj_mat = proj_mat

    def project(self, densities):
        """
        Project the given 3D (volume) densities to 2D (surface) densities,
        using the calculated projection matrix: 'proj_mat'.
        """
        densities = np.array(densities)
        if self.rout.shape != densities.shape:
            raise ValueError("different shapes of rout and given densities")
        return self.proj_mat.dot(densities.T)

    def deproject(self, densities):
        """
        Revert the projection procedure, i.e., deproject the given 2D
        (surface) densities to derive the 3D (volume) densities.

        \curl{N}(R_{i-1}, R_i) = N(R_{i-1}, R_i) * \pi * (R^2_i - R^2_{i-1})

        n(R_{i-1}, R_i) =
            (N(R_{i-1}, R_i) * \pi * (R^2_i - R^2_{i-1}) /
             V_int(R_{i-1}, R_i; R_{i-1}, R_i)) -
            \sum_{j=i+1}^{m} (n(R_{j-1}, R_j) *
                              V_int(R_{j-1}, R_j; R_{i-1}, R_i) /
                              V_int(R_{i-1}, R_i; R_{i-1}, R_i))

        Reference: ref.[2], eq.(A2)
        """
        densities = np.array(densities)
        if self.rout.shape != densities.shape:
            raise ValueError("different shapes of rout and given densities")
        n_3d = np.zeros(densities.shape)
        # peel the onion: from outside inward
        for i in reversed(range(self.N)):
            rin = self.rin[i]
            rout = self.rout[i]
            area = np.pi * (rout**2 - rin**2)
            v_int = self.intersection_volume(rin, rout, rin, rout)
            n_3d[i] = densities[i] * area / v_int
            # subtract the projections from the outer shells
            for j in range(i+1, self.N):
                rin2 = self.rin[j]
                rout2 = self.rout[j]
                v_int2 = self.intersection_volume(rin2, rout2, rin, rout)
                n_3d[i] -= n_3d[j] * v_int2 / v_int
        return n_3d

    @staticmethod
    def intersection_volume(r1, r2, R1, R2):
        """
        Calculate the volume of intersection between the spherical shell of
        r1 <= r <= r2 and the cylinder of R1 <= R <= R2.

        Reference: ref.[2], eq.(A1)
        """
        def trunc_pow(x, p):
            if x <= 0.0:
                return 0
            else:
                return x ** p
            #
        v_int = (4.0*np.pi/3.0) * (trunc_pow((r2**2 - R1**2), 1.5) -
                                   trunc_pow((r2**2 - R2**2), 1.5) +
                                   trunc_pow((r1**2 - R2**2), 1.5) -
                                   trunc_pow((r1**2 - R1**2), 1.5))
        return v_int


def testProjection():
    rout = np.array([1, 2, 3, 4, 5], dtype=float)
    proj = Projection(rout)
    n1 = np.array([1, 1, 1, 1, 1], dtype=float)
    np.testing.assert_array_almost_equal(proj.deproject(proj.project(n1)), n1)
    s2 = np.array([1, 1, 1, 1, 1], dtype=float)
    np.testing.assert_array_almost_equal(proj.project(proj.deproject(s2)), s2)
    print("All tests PASSED!")


class ABModel:
    """
    AB model is a modified version of the beta model, which can roughly
    fit both centrally peaked and cored models, e.g., central excess emission.

    References:
    [1] Pratt & Arnaud, 2002, A&A, 394, 375; eq.(2)
    [2] ref.[1], eq.(10)
    """
    name = "AB model"
    # model parameters
    params = lmfit.Parameters()
    params.add_many(  # (name, value, vary, min, max, expr)
                    ("A",     1.0e-9, True, 0.0, 1.0e-5, None),
                    ("alpha", 0.7,    True, 0.1, 1.1,    None),
                    ("rc",    30.0,   True, 1.0, 1.0e4,  None),
                    ("beta",  0.7,    True, 0.3, 1.1,    None))
    # optimization method
    fit_method = "lbfgsb"
    # whether the 'ydata' and 'yerr' to be scaled in order to reduce
    # the dynamical range for a more stable fitting
    scale = False
    scale_factor = 1.0

    def __init__(self, fit_method="lbfgsb", params=None, scale=True):
        self.fit_method = fit_method
        if params is not None:
            self.load_params(params)
        self.scale = scale

    def load_data(self, xdata, ydata=None, xerr=None, yerr=None,
                  update_params=False):
        self.reset()
        if xdata.ndim == 2 and xdata.shape[1] == 4:
            # 4-column data
            self.xdata = xdata[:, 0].copy()
            self.xerr = xdata[:, 1].copy()
            self.ydata = xdata[:, 2].copy()
            self.yerr = xdata[:, 3].copy()
        else:
            self.xdata = np.array(xdata)
            self.ydata = np.array(ydata)
            self.xerr = np.array(xerr)
            self.yerr = np.array(yerr)
        self.scale_data(update_params=update_params)

    def scale_data(self, update_params=False):
        if self.scale:
            y_min = np.min(self.ydata)
            y_max = np.max(self.ydata)
            self.scale_factor = np.exp(np.log(y_min*y_max) / 2)
            self.ydata /= self.scale_factor
            self.yerr /= self.scale_factor
            if update_params:
                A_min = 1.0
                A_max = np.max(self.ydata)
                self.set_param("A", value=(A_min+A_max)*0.5,
                               min=A_min, max=A_max)

    def reset(self):
        self.fitter = None
        self.fitted = None

    def f_residual(self, params):
        if self.yerr is None:
            return self.model(self.xdata, params) - self.ydata
        else:
            return (self.model(self.xdata, params) - self.ydata) / self.yerr

    def fit(self, method=None):
        if method is None:
            method = self.fit_method
        self.fitter = lmfit.Minimizer(self.f_residual, self.params)
        self.fitted = self.fitter.minimize(method=method)
        self.load_params(self.fitted.params)

    @staticmethod
    def model(x, params):
        parvals = params.valuesdict()
        A = parvals["A"]
        alpha = parvals["alpha"]
        rc = parvals["rc"]
        beta = parvals["beta"]
        return (A * np.power(x/rc, -alpha) *
                np.power((1 + (x/rc)**2), -1.5*beta + 0.5*alpha))

    def f(self, x):
        return self.model(x, self.params) * self.scale_factor

    def get_param(self, name=None):
        """
        Return the requested 'Parameter' object or the whole
        'Parameters' object of no name supplied.
        """
        try:
            return self.params[name]
        except KeyError:
            return self.params

    def set_param(self, name, *args, **kwargs):
        """
        Set the properties of the specified parameter.
        """
        param = self.params[name]
        param.set(*args, **kwargs)

    def dump_params(self, serialize=True):
        """
        Dump the current values/settings for all model parameters,
        and these dumped results can be later loaded by 'load_params()'.
        """
        if serialize:
            return self.params.dumps()
        else:
            return self.params.copy()

    def load_params(self, params):
        """
        Load the provided parameters values/settings.
        """
        if isinstance(params, lmfit.parameter.Parameters):
            self.params = params.copy()
        else:
            p = lmfit.parameter.Parameters()
            p.loads(params)
            self.params = p


class DeprojectSBP:
    """
    Deproject the observed SBP to derive the 3D emission measure (EM)
    profile, using a regularization technique.

    References: ref.[1]
    """
    # input SBP data: [r, r_err, s, s_err]
    r = None
    r_err = None
    s = None
    s_err = None
    # 'Projection' instance for this SBP
    projector = None
    # 'ABModel' instance to fit the deprojected EM profile for rescaling data
    abmodel = None
    # smoothing parameter to balance between fidelity (chisq) and
    # consistency with the applied regularization constraint.
    lbd = 1.0
    # optimization method for scipy minimize
    opt_method = "Powell"
    # scipy optimize results from 'self.deproject()'
    deproject_res = None
    ####
    debug_xlist = []

    def __init__(self, r, r_err=None, s=None, s_err=None,
                 lbd=1.0, opt_method="Powell"):
        self.load_data(r=r, r_err=r_err, s=s, s_err=s_err)
        self.projector = Projection(rout=self.r+self.r_err)
        self.abmodel = ABModel(scale=True)
        self.lbd = lbd
        self.opt_method = opt_method
        self.debug_xlist = []

    def load_data(self, r, r_err=None, s=None, s_err=None):
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

    def deproject(self, lbd=None, opt_method=None):
        """
        Deproject the observed SBP to derive the 3D EM profile by
        minimizing the objective function.
        """
        def fobj(x):
            return self.f_objective(x, scaled=True)

        def callback(x):
            # NOTE: 'x' here is the scaled EM solution
            x_unscaled = self.unscale_data(x)
            self.update_abmodel(x_unscaled)

        if lbd is not None:
            self.lbd = lbd
        if opt_method is None:
            opt_method = self.opt_method
        # initial guess
        em0 = self.projector.deproject(self.s)
        # scale the EM data to reduce the dynamical range
        em0_scaled = self.scale_data(em0, update_params=True)
        # DEBUG
        self.debug_xlist.append(em0_scaled)
        res = scipy.optimize.minimize(fun=fobj, x0=em0_scaled,
                                      method=opt_method,
                                      callback=callback,
                                      options={"disp": True})
        self.deproject_res = res
        self.em = self.unscale_data(res.x)
        return self.em

    def update_abmodel(self, x, xerr=None, update_params=False):
        """
        Load the supplied data into self.abmodel, and perform fitting.

        If the errors/uncertainties is not specified, it is assumed
        to have the same relative errors as the observed SBP.
        """
        self.debug_xlist.append(x)
        if xerr is None:
            x_err = x * self.s_err / self.s
        self.abmodel.load_data(xdata=self.r, xerr=self.r_err,
                               ydata=x, yerr=x_err,
                               update_params=update_params)
        self.abmodel.fit()

    def scale_data(self, x, xerr=None, update_params=False):
        """
        Scale the data (i.e., 3D EM profile) by dividing the fitted
        AB model.

        If the errors/uncertainties is not specified, it is assumed
        to have the same relative errors as the observed SBP.
        """
        self.update_abmodel(x=x, xerr=xerr, update_params=update_params)
        x_fitted = self.abmodel.f(self.r)
        x_scaled = x / x_fitted
        return x_scaled

    def unscale_data(self, x):
        """
        Undo the data scaling by multiplying the same fitted model
        previously used to scale the data.
        """
        x_fitted = self.abmodel.f(self.r)
        x_unscaled = x * x_fitted
        return x_unscaled

    def f_objective(self, x, scaled=False):
        """
        The objective function to be minimized, in order to derive the
        best solution (i.e., deprojected SBP) for the observed SBP.

        This objective function is a combination of plain chi-squared
        and a regularization constraint.

        'lbd' is the parameter to balance the goodness-of-fit and the
        regularization constraint.

        References: ref.[1], eq.(2)
        """
        return (self.f_chisq(x, scaled=scaled) +
                self.lbd * self.f_constraint(x, scaled=scaled))

    def f_residual(self, x, scaled=False):
        """
        Calculate the residuals of each data point for the solution.

        The current solution (i.e., 3D EM profile) is first projected
        into the 2D SBP, then compared to the observed SBP.
        """
        if scaled:
            x = self.unscale_data(x)
        x_2d = self.projector.project(x)
        residuals = (x_2d - self.s) / self.s_err
        return residuals

    def f_chisq(self, x, scaled=False):
        """
        Function to calculate the chi-squared value of the current
        solution with respect to the data.
        """
        chisq = np.sum(self.f_residual(x, scaled=scaled) ** 2)
        return chisq

    def f_constraint(self, x, scaled=False):
        """
        Function to calculate the value of regularization constraint.

        References: ref.[1], eq.(3)
        """
        if not scaled:
            x = self.scale_data(x)
        constraint = np.sum((x[:-1] + x[1:]) ** 2)
        return constraint


def main():
    parser = argparse.ArgumentParser(
            description="Deproject the surface brightness profile (SBP)")
    args = parser.parse_args()


if __name__ == "__main__":
    main()
