# -*- mode: python -*-
#
# Weitian LI
# Created: 2016-06-24
# Updated: 2016-06-24
#

"""
This module contains the parameters/constants used in astronomy
and astrophysics.
"""

import astropy.units as au


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
