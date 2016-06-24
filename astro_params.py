# -*- mode: python -*-
#
# Weitian LI
# Created: 2016-06-24
# Updated: 2016-06-24
#
# Change logs:
# 2016-06-24:
#   * Add class 'ChandraPixel', moved from 'deproject_sbp.py'
#

"""
This module contains the parameters/constants used in astronomy
and astrophysics.
"""

import astropy.units as au
from astropy.cosmology import FlatLambdaCDM


class AstroParams:
    """
    The parameters/constants used in astronomy.

    References:
    [1] ref. [4], eq.(9) below
    """
    # Hubble constant at z=0
    H0 = 71.0  # [ km/s/Mpc ]
    # density of non-relativistic matter in units of the critical density
    # at z=0
    OmegaM0 = 0.27
    # ratio of electron density (n_e) to proton density (n_p) [1]
    ratio_ne_np = 1.211
    # molecular weight per electron (0.3 solar abundance; grsa table) [1]
    mu_e = 1.155
    # atomic mass unit
    m_atom = au.u.to(au.g)  # [ g ]


class ChandraPixel:
    """
    Chandra pixel unit conversions.
    """
    angle = 0.492 * au.arcsec
    z = None
    # cosmology calculator
    cosmo = None
    # angular diameter distance
    D_A = None
    # length of one pixel at the given redshift
    length = None

    def __init__(self, z=None):
        self.z = z
        self.cosmo = FlatLambdaCDM(H0=AstroParams.H0,
                                   Om0=AstroParams.OmegaM0)
        if z is not None:
            self.D_A = self.cosmo.angular_diameter_distance(z)
            self.length = self.D_A * self.angle.to(au.radian).value

    def get_angle(self):
        return self.angle

    def get_length(self, z=None):
        if z is None:
            length = self.length
        else:
            D_A = self.cosmo.angular_diameter_distance(z)
            length = D_A * self.angle.to(au.radian).value
        return length
