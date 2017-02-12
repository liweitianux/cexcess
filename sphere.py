# Copyright (c) 2017 Aaron LI
# MIT license

"""
Spherical utilities.
"""

import numpy as np


def central_angle(p0, points):
    """
    Calculate the central angle(s) between the points ``p0`` with respect to
    the other point(s) ``points`` on the sphere.

    Parameters
    ----------
    p0 : 2-element float tuple/list
        (longitude/R.A., latitude/Dec.) coordinate of the reference point.
        (Unit: deg)
    points : 2-element float tuple/list, or 2-column float `~numpy` array
        Coordinates of the other point(s)
        (Unit: deg)

    Returns
    -------
    angle : float, or 1D float `~numpy` array
        Calculated central angle(s) (Unit: deg)

    Algorithm
    ---------
    (radial, azimuthal, polar): (r, θ, φ)
    central_angle: α
    longitude/R.A.: λ = θ
    latitude/Dec.: δ = π/2 - φ

    Unit vector:
        r1_vec = (cos(θ1)*sin(φ1), sin(θ1)*sin(φ1), cos(φ1))
               = (cos(λ1)*cos(δ1), sin(λ1)*cos(δ1), sin(δ1))
        r2_vec = (cos(θ2)*sin(φ2), sin(θ2)*sin(φ2), cos(φ2))
               = (cos(λ2)*cos(δ2), sin(λ2)*cos(δ2), sin(δ2))

    Therefore the angle (α) between r1_vec and r2_vec:
        cos(α) = r1_vec * r2_vec
               = cos(δ1)*cos(δ2)*cos(λ1-λ2) + sin(δ1)*sin(δ2)

    References
    ----------
    [1] Spherical Coordinates - Wolfram MathWorld
        http://mathworld.wolfram.com/SphericalCoordinates.html
        Eq.(19)
    [2] Great Circle - Wolfram MathWorld
        http://mathworld.wolfram.com/GreatCircle.html
        Eqs.(1,2,4)
    """
    lon0, lat0 = np.deg2rad(p0)
    points = np.deg2rad(points)
    try:
        lon, lat = points[:, 0], points[:, 1]
    except IndexError:
        lon, lat = points  # single point
    prod = (np.cos(lat0) * np.cos(lat) * np.cos(lon0-lon) +
            np.sin(lat0) * np.sin(lat))
    alpha = np.arccos(prod)
    return np.rad2deg(alpha)
