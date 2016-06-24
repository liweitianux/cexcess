#!/usr/bin/env python3
#
# Weitian LI
# Created: 2016-06-10
# Updated: 2016-06-24
#

"""
Project the 3D volume density to 2D surface density and vice versa.

References:
[1] McLaughlin, 1999, ApJ, 117, 2398-2427
"""

import numpy as np


class Projection:
    """
    Class that deals with projection from 3D volume density to 2D
    surface density and vice versa.

    NOTE:
    * The inner-most shell/cylinder is assumed to at the center with inner
      radius of ZERO.
    * Uniform background should be subtracted before carrying out the
      deprojection.
    * The surface density is assumed to be cut at the largest available
      radius, i.e., it is assumed that there isn't any density distributed
      beyond the outer-most shell/cylinder.
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


if __name__ == "__main__":
    testProjection()
