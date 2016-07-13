# -*- mode: python -*-
#
# Aaron LI
# Created: 2016-07-10
# Updated: 2016-07-13
#
# Change logs:
# 2016-07-13:
#   * Improve the `np.array` usage a bit
#

"""
This module implements the spline functions used to interpolate
and/or extrapolate a group of discrete data.

* smoothing spline: R's mgcv::gam()
"""

import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr


class Spline:
    """
    Meta class to fit a spline to the input data.
    """
    # input data
    x = None
    y = None
    # weights of each data point when fitting the spline
    weights = None
    # specifies which axis should be transformed to logarithmic scale
    log10 = None
    # fitted spline for interpolation and extrapolation
    spline = None

    def __init__(self, x, y, weights=None):
        self.x = np.array(x)
        self.y = np.array(y)
        if weights is None:
            self.weights = np.ones(self.x.shape)
        else:
            self.weights = np.array(weights)
        self.log10 = []

    def fit(self, log10=[]):
        """
        The parameter `log10` specifies the axis to be transformed
        to the logarithmic scale before fitted by the spline, e.g.,
        `log10=["x"]`, `log10=["x", "y"]`.
        """
        raise NotImplementedError

    def eval(self, x):
        """
        Evaluate the specified spline at the supplied positions.
        Also check whether the spline was fitted in the log-scale space,
        and transform the evaluated values if necessary.
        """
        raise NotImplementedError


class SmoothSpline(Spline):
    """
    Fit a smoothing spline to the data.

    Currently, the penalized smoothing spline from R's `mgcv::gam()`
    is employed.  In addition, the fitted spline allows extrapolation.
    """
    mgcv = importr("mgcv")

    def fit(self, log10=[]):
        self.log10 = list(map(str.upper, log10))
        if "X" in self.log10:
            x = ro.FloatVector(np.log10(self.x))
        else:
            x = ro.FloatVector(self.x)
        if "Y" in self.log10:
            y = ro.FloatVector(np.log10(self.y))
        else:
            y = ro.FloatVector(self.y)
        weights = ro.FloatVector(self.weights)
        self.spline = self.mgcv.gam(
            ro.Formula("y ~ s(x, bs='ps')"),
            data=ro.DataFrame({"x": x, "y": y}),
            weights=weights, method="REML")

    def eval(self, x):
        x = np.array(x, ndmin=1)
        if "X" in self.log10:
            x_new = ro.ListVector({"x": ro.FloatVector(np.log10(x))})
        else:
            x_new = ro.ListVector({"x": ro.FloatVector(x)})
        y_pred = self.mgcv.predict_gam(self.spline, newdata=x_new)
        if "Y" in self.log10:
            y_pred = 10 ** np.array(y_pred)
        else:
            y_pred = np.array(y_pred)
        #
        if len(y_pred) == 1:
            return y_pred[0]
        else:
            return y_pred
