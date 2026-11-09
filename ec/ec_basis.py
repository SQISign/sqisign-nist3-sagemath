from sage.all import ZZ, matrix, Zmod

import crypto.utils as utils
from ec.ec_point import ECPoint

class ECBasis:
    r"""
    Implements the type ECBasis [Definition 4.4.1].

    A torsion basis `(P, Q)` of `E[2^e]` in x-only coordinates.

    INPUT:

    - ``P``, ``Q`` -- points of type :class:`ECPoint`
    - ``PQ`` -- the point `P - Q`
    """

    def __init__(self, P, Q, PQ):
        if not all(isinstance(R, ECPoint) for R in (P, Q, PQ)):
            raise ValueError("The basis entries must be points")

        self._E = P._parent
        self._P = P
        self._Q = Q
        self._PQ = PQ

    def __repr__(self):
        return f"{self._P},{self._Q},{self._PQ}"

    def generators(self):
        r"""
        Return the list ``[P, Q, P - Q]``.
        """
        return [self._P, self._Q, self._PQ]

    def __iter__(self):
        for _P in self.generators():
            yield _P

    def _biquadratics(self, x1, z1, x2, z2):
        r"""
        [Algorithm 4.91] Biquadratics, on coordinates.
        """
        A = self._E._A
        C = self._E._C

        Bxx = C * (x1 * x2 - z1 * z2) ** 2
        Bxz = C * (x1 * x2 + z1 * z2) * (x1 * z2 + z1 * x2) + 2 * A * x1 * x2 * z1 * z2
        Bzz = C * (x1 * z2 - z1 * x2) ** 2

        return Bxx, Bxz, Bzz

    def _shared_difference(self, x1, z1, x2, z2, x1Q, z1Q, x2Q, z2Q):
        r"""
        [Algorithm 4.93] SharedDifference, on coordinates.
        """

        Bxx, Bxz, Bzz = self._biquadratics(x1, z1, x2, z2)
        Bxx_, Bxz_, Bzz_ = self._biquadratics(x1Q, z1Q, x2Q, z2Q)

        x = 2 * (Bxz_ * Bxx - Bxz * Bxx_)
        z = Bzz_ * Bxx - Bzz * Bxx_

        return x, z

    def _diff_points(self, basis2):
        r"""
        [Algorithm 4.105] DiffPoints.

        INPUT:

        - ``basis2`` -- an :class:`ECBasis` ``(R, S, R - S)``

        OUTPUT:

        A tuple ``(R - P, S - P, R - Q, S - Q)``.
        """

        xP, zP = self._P.xz()
        xQ, zQ = self._Q.xz()
        xPQ, zPQ = self._PQ.xz()

        R, S, RS = basis2
        xR, zR = R.xz()
        xS, zS = S.xz()
        xRS, zRS = RS.xz()

        Bxx, Bxz, Bzz = self._biquadratics(xR, zR, xP, zP)

        xRP = Bxz + utils.fp2_sqrt(Bxz**2 - Bxx * Bzz)
        zRP = Bzz

        xSP, zSP = self._shared_difference(xS, zS, xP, zP, xRP, zRP, xRS, zRS)
        xRQ, zRQ = self._shared_difference(xR, zR, xQ, zQ, xRP, zRP, xPQ, zPQ)
        xSQ, zSQ = self._shared_difference(xS, zS, xQ, zQ, xSP, zSP, xPQ, zPQ)

        RP = ECPoint(self._E, [xRP, zRP])
        SP = ECPoint(self._E, [xSP, zSP])
        RQ = ECPoint(self._E, [xRQ, zRQ])
        SQ = ECPoint(self._E, [xSQ, zSQ])

        return RP, SP, RQ, SQ

    def change_of_basis(self, basis2, e):
        r"""
        [Algorithm 4.106] ChangeOfBasis.

        INPUT:

        - ``basis2`` -- an :class:`ECBasis` of `E[2^e]`
        - ``e`` -- an integer with `e \le f`

        OUTPUT:

        A `2 \times 2` matrix ``M`` over `\ZZ/2^e\ZZ` such that
        ``basis2`` `= B \cdot 2^{f - e} M`, with `B` this basis of `E[2^f]`.
        """
        f = (self._E._p + 1).valuation(2)

        P = self._P
        Q = self._Q
        PQ = self._PQ

        R, S, _ = basis2

        RP, SP, RQ, SQ = self._diff_points(basis2)

        # no need for normalisation, Tate pairing already does it for you
        zeta = P.tate_pairing(Q, PQ, f) ** (2 ** (f-e))
        zeta_i = [
            R.tate_pairing(Q, RQ, e),
            S.tate_pairing(Q, SQ, e),
            1 / R.tate_pairing(P, RP, e),
            1 / S.tate_pairing(P, SP, e)
        ]

        x_i = [utils.normalized_dlog(z_i, 1/zeta, e) for z_i in zeta_i]
        return matrix(Zmod(2**e), 2, x_i)

    def act_on_basis(self, M):
        r"""
        [Algorithm 3.7] ActOnBasis.

        INPUT:

        - ``M`` -- a `2 \times 2` matrix `\begin{pmatrix} x_1 & x_2 \\ x_3 & x_4
          \end{pmatrix}` over `\ZZ/2^e\ZZ`

        OUTPUT:

        The :class:`ECBasis` `([x_1]P + [x_3]Q, [x_2]P + [x_4]Q)`.
        """
        m1, m3, m2, m4 = map(ZZ, M.list())
        mod = M.base_ring().characteristic()
        P, Q, PQ = self
        P_eval = P.ec_ladder_biscalar(Q, PQ, m1, m2)
        Q_eval = P.ec_ladder_biscalar(Q, PQ, m3, m4)
        PQ_eval = P.ec_ladder_biscalar(Q, PQ, (m1-m3) % mod, (m2-m4)%mod)

        return ECBasis(P_eval, Q_eval, PQ_eval)
