from sage.all import EllipticCurve, kronecker

import crypto.utils as utils
from ec.ec_point import ECPoint
from ec.ec_basis import ECBasis

class ECCurve:
    r"""
    Implements the type ECCurve [Section 3.1.2].

    The Montgomery curve `y^2 = x^3 + Ax^2 + x` over `\mathbb{F}_{p^2}`.

    INPUT:

    - ``A`` -- an element of `\mathbb{F}_{p^2}` with `A^2 \neq 4`
    """

    def __init__(self, A):
        self._A = A
        self._Fp2 = self._A.parent()
        self._p = self._Fp2.characteristic()

        if self._A**2 == self._Fp2(4):
            raise ValueError("The coefficients do not describe a Montgomery curve")

        self._i = self._Fp2.gen()
        if self._i**2 != self._Fp2(-1):
            raise ValueError(
                "The base field is not a quadratic extension of a prime field Fp, where p = 3 mod 4"
            )

        self.E = EllipticCurve([0, self._A, 0, 1, 0])

        self._C = self._Fp2(1)
        self._A24 = self._A + 2 * self._C
        self._C24 = 4 * self._C

    def __eq__(self, other):
        if not isinstance(other, ECCurve):
            return False
        return self.A_coefficient() == other.A_coefficient()

    def A_coefficient(self):
        return self._A

    def j_inv(self):
        return self.E.j_invariant()

    def random_point(self):
        r"""
        Return a random point on the curve, in `(x : 1)` coordinates.
        """
        A = self._A
        x = self._Fp2.random_element()
        y_sqrd = x**3 + A * x**2 + x
        while not y_sqrd.is_square():
            x = self._Fp2.random_element()
            y_sqrd = x**3 + A * x**2 + x

        return ECPoint(self, [x, 1])

    def torsion_basis_to_hint(self, e, eh):
        r"""
        [Algorithm 4.94] TorsionBasisToHint.

        Compute a deterministic basis of `E[2^e]` and the hint to recover it.

        INPUT:

        - ``e`` -- an integer with `e \le f`
        - ``eh`` -- an integer with `eh < e`

        OUTPUT:

        A list ``[B, (h, hA)]`` with ``B`` an :class:`ECBasis` of `E[2^e]` and
        ``(h, hA)`` the hint, with ``h`` set to `0` if it exceeds `127`.
        """

        A = self._A
        p = self._p
        i = self._i
        two_e = 2**e
        cofactor = (p + 1) // two_e

        if A.is_square():
            hA = 1
        else:
            hA = 0

        h = 0
        if hA == 1:
            while True:
                h += 1
                xR = -A
                zR = 1 + i * h
                x = xR / zR
                x_sqrd = x**2
                y_sqrd = x_sqrd * x + A * x_sqrd + x
                if kronecker(h**2 + 1, p) == -1 and y_sqrd.is_square():
                    break
        else:
            while True:
                h += 1
                xR = h * A  # we should use addition instead, just a microsaving
                zR = 1
                y_sqrd = xR**3 + A * xR**2 + xR
                if y_sqrd.is_square():
                    break

        S = ECPoint(self, [-(xR + A * zR), zR]).ec_ladder(cofactor)
        R = ECPoint(self, [xR, zR]).ec_ladder(cofactor)
        RS = R.ec_projective_difference(S)

        if eh > 0:
            two_e_eh = 2 ** (e - eh)
            RS2 = R.ec_x_add(S, RS)
            T1 = RS.ec_ladder(two_e_eh)
            T2 = RS2.ec_ladder(two_e_eh)

            X_T1, Z_T1 = T1.xz()
            X_T2, Z_T2 = T2.xz()
            Z = (X_T1 * Z_T2 - Z_T1 * X_T2) * (Z_T1 * Z_T2).conjugate() * Z_T1 * (Z_T1.conjugate()) * Z_T2 * (Z_T2.conjugate())
            if utils.fp2_sqrt(Z**2) == Z:
                RS = RS2

        if h > 127:
            h = 0

        S, RS = RS, S
        return [ECBasis(R, S, RS), (h, hA)]

    def torsion_basis_from_hint(self, e, hint):
        r"""
        [Algorithm 4.95] TorsionBasisFromHint.

        INPUT:

        - ``e`` -- an integer with `e \le f`
        - ``hint`` -- a hint ``(h, hA)`` as returned by
          :meth:`torsion_basis_to_hint`

        OUTPUT:

        The :class:`ECBasis` of `E[2^e]` associated to ``hint``.
        """

        A = self._A
        p = self._p
        i = self._i
        two_e = 2**e
        cofactor = (p + 1) // two_e

        h, hA = hint

        if h == 0:
            basis, _ = self.torsion_basis_to_hint(self, e)
            return basis

        if hA == 1:
            xR = -A
            zR = 1 + i * h
        else:
            xR = h * A
            zR = 1

        S = ECPoint(self, [-(xR + A * zR), zR]).ec_ladder(cofactor)
        R = ECPoint(self, [xR, zR]).ec_ladder(cofactor)
        RS = R.ec_projective_difference(S)
        S, RS = RS, S

        return ECBasis(R, S, RS)

    def __call__(self, P_coords):
        r"""
        Return the point with coordinates ``P_coords = [X, Z]``.
        """
        return ECPoint(self, P_coords)

    def __repr__(self):
        return self.E.__repr__()

    def zero(self):
        r"""
        Return the point at infinity `(1 : 0)`.
        """
        return ECPoint(self, [self._Fp2.one(), self._Fp2.zero()])
