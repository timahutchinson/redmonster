#
# grid_spline.py
#
# Code for one-dimensional cubic splines on a
# uniform grid, including analytic slope, curvature,
# and extremum evaluation.
#
# Most convenient interface is via GridSpline class,
# which encapsulates the low-level routines.
#
# A. Bolton, U. of Utah, 2010-2014
#
from __future__ import division

import numpy as n

def tri_diag(a, b, c, r):
    """
    Tri-diagonal solver
    """
    ndim = len(b)
    alpha = a.copy()
    beta = b.copy()
    for i in range(1,ndim):
        beta[i] = b[i] - alpha[i] * c[i-1] / beta[i-1]
    gamma = c / beta
    y = r.copy()
    y[0] = r[0] / beta[0]
    for i in range(1,ndim):
        y[i] = (r[i] - alpha[i] * y[i-1]) / beta[i]
    x = y.copy()
    x[ndim-1] = y[ndim-1]
    for i in range(ndim-2, -1, -1):
        x[i] = y[i] - gamma[i] * x[i+1]
    return x

def spline_get_ms(y):
    """
    Compute knot slopes to initialize spline
    """
    bign = len(y) - 1
    m = 0. * y
    m[0] = 2.0 * y[1] - 1.5 * y[0] - 0.5 * y[2]
    m[bign] = 1.5 * y[bign] + 0.5 * y[bign-2] - 2.0 * y[bign-1]
    r = 3.0 * (y[2:] - y[:-2])
    r[0] = r[0] - m[0]
    r[bign-2] = r[bign-2] - m[bign]
    on_diag = n.zeros(bign-1) + 4.0
    off_diag = n.zeros(bign-1) + 1.0
    m[1:bign] = tri_diag(off_diag, on_diag, off_diag, r)
    return m

def spline_get_val(y, m, x):
    """
    Evaluate spline value at positions x
    """
    intervals = len(y) - 1
    i = n.int32(x) + 1
    # (the following is a hack to keep the upper
    # bound in the valid interval range:)
    i = i - i // (intervals + 1)
    d1 = x - i + 1.
    d2 = d1 - 1.0
    spline_val = (y[i] * d1**3 - y[i-1] * d2**3
                  + (m[i] - 3.0 * y[i]) * d1**2 * d2
                  + (m[i-1] + 3.0 * y[i-1]) * d1 * d2**2)
    return spline_val

def spline_get_slope(y, m, x):
    """
    Evaluate spline slope at positions x
    """
    intervals = len(y) - 1
    i = n.int32(x) + 1
    # (the following is a hack to keep the upper
    # bound in the valid interval range:)
    i = i - i // (intervals + 1)
    d1 = x - i + 1.
    d2 = d1 - 1.0
    spline_slope = (m[i] * d1**2 + m[i-1] * d2**2 +
                    (2.0*m[i] + 2.0*m[i-1] - 6.0*y[i] + 6.0*y[i-1]) * d1 * d2)
    return spline_slope

def spline_get_curv(y, m, x):
    """
    Evaluate spline curvature at positions x
    """
    intervals = len(y) - 1
    i = n.int32(x) + 1
    # (the following is a hack to keep the upper
    # bound in the valid interval range:)
    i = i - i // (intervals + 1)
    d1 = x - i + 1.
    d2 = d1 - 1.0
    spline_curv = (2.0 * m[i] * d1 + 2.0 * m[i-1] * d2 +
                   (2.0*m[i] + 2.0*m[i-1] - 6.0*y[i] + 6.0*y[i-1]) * (d1 + d2))
    return spline_curv

def spline_get_max(y, m):
    """
    Find positions of analytic maxima of spline
    """
    bign = len(y) - 1
    xval = n.zeros(bign) - 1.0
    #Quadratic derivative coefficients in the intervals:
    a = (3.0 * (m[0:bign] + (n.roll(m, -1))[0:bign])
         + 6.0 * (y[0:bign] - (n.roll(y, -1))[0:bign]))
    b = (-2.0 * (2.0 * m[0:bign] + (n.roll(m, -1))[0:bign]
                 + 3.0 * (y[0:bign] - (n.roll(y, -1))[0:bign])))
    c = (m[0:bign])
    # Discriminant:
    d = b**2 - 4.0 * a * c
    # Find any linear-root maxima:
    lroots = (n.where((a == 0) * (b < 0)))[0]
    if len(lroots) > 0:
        xval[lroots] = -c[lroots] / b[lroots]
    # Find any quadratic-root maxima:
    qroots = (n.where((a != 0) * (d > 0)))[0]
    if len(qroots) > 0:
        xval[qroots] = -0.5 * (b[qroots] + n.sqrt(d[qroots])) / a[qroots]
    # Find roots that are within the necessary interval bounds:
    roots = (n.where((xval >= 0.0) * (xval < 1.0)))[0]
    if len(roots) <= 0:
        return n.asarray([])
    # Transform root values to global x-coordinate and return.
    xval = xval + n.arange(bign)
    xval = xval[roots]
    return xval

# OOP interface to this business:
class GridSpline:
    """
    Initialize a spline object for uniformly gridded 1D data.

    Calling syntax:
      GS = GridSpline(y)

    where y is an array of values to be splined.

    The abscissa for the spline is taken to be a zero-based
    vector of integers of length equal to the y-vector.

    A. Bolton, U. of Utah, 2010-2014
    """
    def __init__(self, y):
        self.y = y.copy()
        self.ms = spline_get_ms(self.y)

    def get_val(self, x):
        """
        Return spline evaluated at abscissa positions x.
        """
        return spline_get_val(self.y, self.ms, x)

    def get_slope(self, x):
        """
        Return analytic derivative of spline evaluated at
        abscissa positions x.
        """
        return spline_get_slope(self.y, self.ms, x)

    def get_curv(self, x):
        """
        Return analytic curvature of spline evaluated at
        abscissa positions x.
        """
        return spline_get_curv(self.y, self.ms, x)

    def get_max(self):
        """
        Return analytically determined locations of maxima
        of spline over domain of original values.
        """
        return spline_get_max(self.y, self.ms)

    def get_min(self):
        """
        Return analytically determined locations of minima
        of spline over domain of original values.
        """
        return spline_get_max(-self.y, -self.ms)
