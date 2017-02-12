# Copyright (c) 2016 Aaron LI
# MIT license
#
# Created: 2016-06-26
#
# Change logs:
# 2016-07-04:
#   * Add "report()" method to class "FittingModel"
#

from collections import OrderedDict

import numpy as np
import lmfit


class FittingModel:
    """
    Base/Meta class for model fitting, with data and parameters scaling.
    """
    name = ""
    params = lmfit.Parameters()
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

    @staticmethod
    def model(x, params):
        pass

    def f(self, x):
        return self.model(x, self.params) * self.scale_factor

    def load_data(self, xdata, ydata=None, xerr=None, yerr=None,
                  update_params=False):
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
        """
        Scale the ydata and yerr to reduce their dynamical ranges,
        for a more stable model fitting.
        """
        if self.scale:
            y_min = np.min(self.ydata)
            y_max = np.max(self.ydata)
            self.scale_factor = np.exp(np.log(y_min*y_max) / 2)
            self.ydata /= self.scale_factor
            self.yerr /= self.scale_factor
            if update_params:
                self.scale_params()

    def scale_params(self):
        """
        Scale the paramters' min/max values accordingly.
        """
        pass

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

    def report(self, rtype):
        """
        Report the fitting results, e.g., g.o.f, chisqr, parameters, etc.
        """
        if rtype == "fitting":
            fitted = self.fitted
            results = OrderedDict([
                ("nfev",   fitted.nfev),
                ("ndata",  fitted.ndata),
                ("nvarys", fitted.nvarys),  # number of variable parameters
                ("nfree",  fitted.nfree),  # degree of freedom
                ("chisqr", fitted.chisqr),
                ("redchi", fitted.redchi),
                ("aic",    fitted.aic),
                ("bic",    fitted.bic),
            ])
        elif rtype == "parameters":
            results = OrderedDict([
                (pn, [par.value, par.min, par.max, par.vary])
                for pn, par in self.params.items()
            ])
        else:
            raise ValueError("invalid rtype: %s" % rtype)
        return results


class ABModel(FittingModel):
    """
    AB model is a modified beta model, which can roughly fit both
    centrally peaked and cored models, e.g., central excess emission.

    This model is used here to constrain the deprojected 3D gas density
    profile, in order to require it is smooth enough.

    References:
    [1] Pratt & Arnaud, 2002, A&A, 394, 375; eq.(2)
    [2] Croston et al. 2006, A&A, 459, 1007-1019; eq.(10)
    """
    name = "AB model"
    # model parameters
    params = lmfit.Parameters()
    params.add_many(  # (name, value, vary, min, max, expr)
                    ("A",     1.0e-9, True, 0.0, 1.0e-5, None),
                    ("alpha", 0.7,    True, 0.1, 1.1,    None),
                    ("rc",    30.0,   True, 1.0, 1.0e4,  None),
                    ("beta",  0.7,    True, 0.3, 1.1,    None))

    def scale_params(self):
        A_min = 1.0
        A_max = np.max(self.ydata)
        self.set_param("A", value=(A_min+A_max)*0.5,
                       min=A_min, max=A_max)

    @staticmethod
    def model(x, params):
        parvals = params.valuesdict()
        A = parvals["A"]
        alpha = parvals["alpha"]
        rc = parvals["rc"]
        beta = parvals["beta"]
        return (A * np.power(x/rc, -alpha) *
                np.power((1 + (x/rc)**2), -1.5*beta + 0.5*alpha))


class PLCModel(FittingModel):
    """
    PLC model consists of a powerlaw and an constant, that is used
    to fit/approximate the outer SBP.
    Therefore, the fitted constant is used to subtract the uniform
    background from the SBP, and the fitted powerlaw index is used
    to extrapolate the SBP in order to mitigate the deprojection
    errors due to FoV limit.

    NOTE:
    I think the uniform background (i.e., by fitting the whole or
    core-excluded SBP) should be subtracted from the SBP first, then
    adopt this PLCModel to fit the outer part of SBP, with the 'bkg'
    parameter fixed at zero.
    """
    name = "PLC model"
    # model parameters
    params = lmfit.Parameters()
    params.add_many(  # (name, value, vary, min, max, expr)
                    ("A",     1.0e-9, True,  0.0, 1.0e-5, None),
                    ("rmin",  30.0,   False, 1.0, 1.0e4,  None),
                    ("alpha", 1.6,    True,  0.4, 2.8,    None),
                    ("bkg",   0.0,    False, 0.0, 1.0e-5, None))

    def load_data(self, xdata, ydata=None, xerr=None, yerr=None,
                  update_params=False):
        super().load_data(xdata=xdata, ydata=ydata, xerr=xerr, yerr=yerr,
                          update_params=update_params)
        self.set_param("rmin", value=np.min(xdata), vary=False)

    def scale_params(self):
        ymin = np.min(self.ydata)
        ymax = np.max(self.ydata)
        self.set_param("A", value=ymax, min=ymax/10.0, max=ymax*10.0)
        self.set_param("bkg", value=ymin, min=0.0, max=ymin)

    @staticmethod
    def model(x, params):
        parvals = params.valuesdict()
        A = parvals["A"]
        rmin = parvals["rmin"]
        alpha = parvals["alpha"]
        bkg = parvals["bkg"]
        return A * np.power(x/rmin, -alpha) + bkg
