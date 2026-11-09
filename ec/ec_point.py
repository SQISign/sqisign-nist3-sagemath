from sage.all import ZZ

import crypto.utils as utils

class ECPoint:
    r"""
    Implements the type ECPoint [Section 3.1.2].

    A point on a Montgomery curve in x-only projective coordinates `(X : Z)`.

    INPUT:

    - ``parent`` -- the :class:`ECCurve` of the point
    - ``coords`` -- the list ``[X, Z]``
    """

    def __init__(self, parent, coords):
        self._x = parent._Fp2(coords[0])
        self._z = parent._Fp2(coords[1])
        self._parent = parent

    def __repr__(self):
        return f"({self._x} : {self._z})"

    def xz(self):
        return self._x, self._z

    def __eq__(self, Q):
        Fp2 = self._parent._Fp2
        zero = Fp2(0)

        x0, z0 = self.xz()
        x1, z1 = Q.xz()

        if z0 == zero:
            return z1 == zero

        if z1 == zero:
            return z0 == zero

        return x0 * z1 == x1 * z0

    def _ec_dbl(self, xP, zP):
        r"""
        [Algorithm 4.84] xDBL, on coordinates.
        """
        Fp2 = self._parent._Fp2
        zero = Fp2(0)

        # identity
        if zP == zero:
            return [1, 0]

        A24 = self._parent._A24
        C24 = self._parent._C24

        t0 = xP + zP
        t0 = t0**2
        t1 = xP - zP
        t1 = t1**2
        t2 = t0 - t1
        t1 = t1 * C24
        x2P = t0 * t1
        t0 = t2 * A24
        t0 = t0 + t1
        z2P = t0 * t2

        return x2P, z2P

    def ec_dbl(self):
        r"""
        [Algorithm 4.84] xDBL.

        OUTPUT:

        The point `[2]P`.
        """

        xP = self._x
        zP = self._z
        x2P, z2P = self._ec_dbl(xP, zP)

        return self._parent([x2P, z2P])

    def _ec_x_add(self, xP, zP, xQ, zQ, xPQ, zPQ):
        r"""
        [Algorithm 4.85] xADD, on coordinates.
        """
        Fp2 = self._parent._Fp2
        zero = Fp2(0)

        # case where P == Q
        if zPQ == zero:
            x, z = self._ec_dbl(xP, zP)
            return x, z

        t0 = xP + zP
        t1 = xP - zP
        t2 = xQ + zQ
        t3 = xQ - zQ
        t0 = t0 * t3
        t1 = t1 * t2
        t2 = t0 + t1
        t3 = t0 - t1
        t2 = t2**2
        t3 = t3**2
        x = zPQ * t2
        z = xPQ * t3

        return x, z

    def ec_x_add(self, Q, PQ):
        r"""
        [Algorithm 4.85] xADD.

        INPUT:

        - ``Q`` -- a point
        - ``PQ`` -- the point `P - Q`

        OUTPUT:

        The point `P + Q`.
        """

        xP = self._x
        zP = self._z

        xQ, zQ = Q.xz()
        xPQ, zPQ = PQ.xz()

        x, z = self._ec_x_add(xP, zP, xQ, zQ, xPQ, zPQ)

        return self._parent([x, z])

    def ec_dbl_iter(self, n):
        r"""
        Iterate :meth:`ec_dbl` ``n`` times.

        OUTPUT:

        The point `[2^n]P`.
        """

        xP = self._x
        zP = self._z

        for _ in range(0, n):
            xP, zP = self._ec_dbl(xP, zP)

        return self._parent([xP, zP])

    def _ec_x_dbl_add(self, xP, zP, xQ, zQ, xPQ, zPQ):
        r"""
        [Algorithm 4.86] xDBLADD, on coordinates.
        """

        Fp2 = self._parent._Fp2
        zero = Fp2(0)

        # case where P == Q
        if zPQ == zero:
            x2P, z2P = self._ec_dbl(xP, zP)
            return x2P, z2P, x2P, z2P

        A24 = self._parent._A24
        C24 = self._parent._C24

        t0 = xP + zP
        t1 = xP - zP
        x2P = t0 ** 2
        t2 = xQ - zQ
        xPQ_plus = xQ + zQ
        t0 = t0 * t2
        z2P = t1 ** 2
        t1 = t1 * xPQ_plus
        t2 = x2P - z2P
        z2P = z2P * C24
        x2P = x2P * z2P
        xPQ_plus = A24 * t2
        zPQ_plus = t0 - t1
        z2P = z2P + xPQ_plus
        xPQ_plus = t0 + t1
        z2P = z2P * t2
        zPQ_plus = zPQ_plus**2
        xPQ_plus = xPQ_plus**2
        zPQ_plus = zPQ_plus * xPQ
        xPQ_plus = xPQ_plus * zPQ

        return x2P, z2P, xPQ_plus, zPQ_plus

    def ec_x_dbl_add(self, Q, PQ):
        r"""
        [Algorithm 4.86] xDBLADD.

        INPUT:

        - ``Q`` -- a point
        - ``PQ`` -- the point `P - Q`

        OUTPUT:

        A tuple ``([2]P, P + Q)``.
        """
        xP = self._x
        zP = self._z

        xQ, zQ = Q.xz()
        xPQ, zPQ = PQ.xz()

        x2P, z2P, xPQ_plus, zPQ_plus = self._ec_x_dbl_add(xP, zP, xQ, zQ, xPQ, zPQ)

        two_P = self._parent([x2P, z2P])
        P_plus_Q = self._parent([xPQ_plus, zPQ_plus])

        return two_P, P_plus_Q

    def _ec_ladder(self, r):
        r"""
        [Algorithm 4.87] Ladder, on coordinates.
        """

        xP = self._x
        zP = self._z

        x0 = self._parent._Fp2(1)
        z0 = self._parent._Fp2(0)
        x1 = xP
        z1 = zP

        binary_r = format(r, "b")
        for k in binary_r:
            if k == "1":
                x1, z1, x0, z0 = self._ec_x_dbl_add(x1, z1, x0, z0, xP, zP)
            else:
                x0, z0, x1, z1 = self._ec_x_dbl_add(x0, z0, x1, z1, xP, zP)

        return x0, z0

    def ec_ladder(self, r):
        r"""
        [Algorithm 4.87] Ladder.

        INPUT:

        - ``r`` -- a non-negative integer

        OUTPUT:

        The point `[r]P`.
        """
        x, z = self._ec_ladder(r)
        return self._parent([x, z])

    def _ec_ladder3pt(self, xP, zP, xQ, zQ, xPQ, zPQ, r):
        r"""
        [Algorithm 4.88] Ladder3pt, on coordinates.
        """

        x0, z0 = xP, zP
        x1, z1 = xQ, zQ
        x2, z2 = xPQ, zPQ

        binary_r = format(r, "b")
        for k in binary_r[::-1]:
            if k == "1":
                x1, z1, x0, z0 = self._ec_x_dbl_add(x1, z1, x0, z0, x2, z2)
            else:
                x1, z1, x2, z2 = self._ec_x_dbl_add(x1, z1, x2, z2, x0, z0)

        return x0, z0

    def ec_ladder3pt(self, Q, PQ, r):
        r"""
        [Algorithm 4.88] Ladder3pt.

        INPUT:

        - ``Q`` -- a point
        - ``PQ`` -- the point `P - Q`
        - ``r`` -- a non-negative integer

        OUTPUT:

        The point `P + [r]Q`.
        """

        xP, zP = self._x, self._z
        xQ, zQ = Q.xz()
        xPQ, zPQ = PQ.xz()

        x, z = self._ec_ladder3pt(xP, zP, xQ, zQ, xPQ, zPQ, r)

        return self._parent([x, z])

    def _ec_ladder_biscalar(self, xP, zP, xQ, zQ, xPQ, zPQ, m, n):
            r"""
            [Algorithm 4.89] LadderBiscalar, on coordinates.
            """

            if m%2 == 0 and n%2 == 1:
                σ0, σ1 = 1, 0
            else:
                σ0, σ1 = 0, 1

            m_, n_ = m , n
            if m%2 == 0:
                m_ = m_ - 1
            if n%2 == 0:
                n_ = n_ - 1

            b0 = format(m_, 'b')
            b1 = format(n_, 'b')

            #computing k
            k_m, k_n = len(b0), len(b1)
            if k_m > k_n:
                b1 = ((k_m - k_n) * "0")+ b1
                k = k_m

            elif k_n > k_m:
                b0 = ((k_n - k_m) * "0")+ b0
                k = k_n

            else:
                k = k_n

            b0 = b0[::-1] + "0"
            b1 = b1[::-1] + "0"

            r = (2*k) * [0]
            for i in range(0, k):

                b0_, b1_ = b0, b1

                if σ0 == 1:
                    b0_ = b1
                    b1_ = b0

                r[2*i] = ZZ(b0_[i]) ^ ZZ(b0_[i+1])
                r[2*i+1] = ZZ(b1_[i]) ^ ZZ(b1_[i+1])

                if r[2*i+1] == 1:
                    t = σ0
                    σ0 = σ1
                    σ1 = t

            xR0, zR0 = 1, 0
            xT0, zT0 = xP, zP
            xT1, zT1 = xQ, zQ

            if σ0 == 0:
                xR1, zR1 = xT0, zT0
                xR2, zR2 = xT1, zT1
            else:
                xR1, zR1 = xT1, zT1
                xR2, zR2 = xT0, zT0

            xD1, zD1 = xR1, zR1
            xD2, zD2 = xR2, zR2

            xR2, zR2 = self._ec_x_add(xR1, zR1, xR2, zR2, xPQ, zPQ)
            xF1, zF1 = xR2, zR2
            xF2, zF2 = xPQ, zPQ

            for i in range(k-1, -1, -1):
                h = r[2*i] + r[2*i+1]
                if h%2 == 0:
                    xT0, zT0 = xR0, zR0
                else:
                    xT0, zT0 = xR1, zR1

                xT0, zT0 = xT0, zT0
                xT1, zT1 = xR2, zR2

                if h < 2:
                    xT0, zT0 = self._ec_dbl(xT0, zT0)
                else:
                    xT0, zT0 = self._ec_dbl(xT1, zT1)

                if r[2*i+1] == 0:
                    xT1, zT1 = xR0, zR0
                    xT2, zT2 = xR1, zR1
                else:
                    xT1, zT1 = xR1, zR1
                    xT2, zT2 = xR2, zR2

                if r[2*i+1] == 1:
                    xTMP, zTMP = xD1, zD1
                    xD1, zD1 = xD2, zD2
                    xD2, zD2 = xTMP, zTMP

                xT1, zT1 = self._ec_x_add(xT1, zT1, xT2, zT2, xD1, zD1)
                xT2, zT2 = self._ec_x_add(xR0, zR0, xR2, zR2, xF1, zF1)

                if h%2 == 1:
                    xTMP, zTMP = xF1, zF1
                    xF1, zF1 = xF2, zF2
                    xF2, zF2 = xTMP, zTMP

                xR0, zR0 = xT0, zT0
                xR1, zR1 = xT1, zT1
                xR2, zR2 = xT2, zT2

            t = (m%2) + (n%2)

            if t == 0:
                x, z = xR0, zR0
            elif t == 1:
                x, z = xR1, zR1
            else:
                x, z = xR2, zR2

            return x, z


    def ec_ladder_biscalar(self, Q, PQ, k_P, k_Q):
        r"""
        [Algorithm 4.89] LadderBiscalar.

        INPUT:

        - ``Q`` -- a point
        - ``PQ`` -- the point `P - Q`
        - ``k_P``, ``k_Q`` -- non-negative integers

        OUTPUT:

        The point `[k_P]P + [k_Q]Q`.
        """

        xP, zP = self._x, self._z
        xQ, zQ = Q.xz()
        xPQ, zPQ = PQ.xz()

        x, z = self._ec_ladder_biscalar(xP, zP, xQ, zQ, xPQ, zPQ, k_P, k_Q)

        return self._parent([x, z])

    def ec_projective_difference(self, Q):
        r"""
        [Algorithm 4.92] ProjectiveDifference.

        INPUT:

        - ``Q`` -- a point

        OUTPUT:

        A deterministic choice of either `P - Q` or `P + Q`.
        """

        xP, zP = self._x, self._z
        xQ, zQ = Q.xz()

        A = self._parent._A
        C = self._parent._C

        Bxx = C * (xP * xQ - zP * zQ) ** 2
        Bxz = C * (xP * xQ + zP * zQ) * (xP * zQ + zP * xQ) + 2 * A * xP * xQ * zP * zQ
        Bzz = C * (xP * zQ - zP * xQ) ** 2
        gamma = C * (C.conjugate() * (xP * zQ - zP * xQ).conjugate())**2

        Bxx = gamma * Bxx
        Bxz = gamma * Bxz
        Bzz = gamma * Bzz
        delta = utils.fp2_sqrt(Bxz**2 - Bxx * Bzz)
        xPQ = Bxz - delta
        zPQ = Bzz
        return self._parent([xPQ, zPQ])

    def ec_bary_components(self, Q, PQ):
        r"""
        [Algorithm 4.107] BARYComponents.

        INPUT:

        - ``Q`` -- a point
        - ``PQ`` -- the point `P - Q`

        OUTPUT:

        A list ``[u, v, w]`` such that `P - Q = (u + v : w)` and
        `P + Q = (u - v : w)`.
        """
        xP, zP = self._x, self._z
        xQ, zQ = Q.xz()
        xPQ, zPQ = PQ.xz()

        t0 = xP + zP
        t1 = xP - zP
        t2 = xQ + zQ
        t3 = xQ - zQ
        t0 = t0 * t3
        t1 = t1 * t2
        t2 = t1 + t0
        t3 = t1 - t0
        t0 = t2 * zPQ
        t0 = t0**2
        t1 = t3 * xPQ
        t2 = t3 * zPQ
        t3 = t1**2
        u = t3 + t0
        v = t3 - t0
        w = t1 * t2
        w = w + w

        return [u, v, w]

    def ec_gluing_ladder(self, Q, PQ, e):
        r"""
        [Algorithm 4.90] GluingLadder.

        INPUT:

        - ``Q`` -- a point
        - ``PQ`` -- the point `P - Q`
        - ``e`` -- a positive integer

        OUTPUT:

        A tuple ``([2^e]P, [2^e - 1]P + Q)``.
        """
        R0 = Q
        R1 = self
        for _ in range(e):
            R1, R0 = R1.ec_x_dbl_add(R0, PQ)
        return (R1, R0)

    @staticmethod
    def _cubical_diff_add(xP, zP, xQ, zQ, inv_xPQ):
        r"""
        [Algorithm 4.100] CubicalDiffAdd, on coordinates.
        """
        a = xP + zP
        b = xP - zP
        c = xQ + zQ
        d = xQ - zQ
        xPQ = (a * d + b * c) ** 2
        zPQ = (a * d - b * c) ** 2
        xPQ = xPQ * inv_xPQ
        return xPQ, zPQ

    def _cubical_dbl(self, xP, zP):
        r"""
        [Algorithm 4.101] CubicalDbl, on coordinates.
        """
        A24 = self._parent._A24
        C24 = self._parent._C24

        a = (xP + zP) ** 2
        b = (xP - zP) ** 2
        c = a - b
        x2 = a * b
        z2 = c * (b + A24 / C24 * c)
        return x2, z2

    def _cubical_ladder(self, e, xPQ, zPQ, xP, zP, xQ):
        r"""
        [Algorithm 4.102] CubicalLadder, on coordinates.
        """
        xnPQ, znPQ = xPQ, zPQ
        xnP, znP = xP, zP
        inv_xQ = xQ ** (-1)

        for _ in range(0, e):
            xnPQ, znPQ = self._cubical_diff_add(xnPQ, znPQ, xnP, znP, inv_xQ)
            xnP, znP = self._cubical_dbl(xnP, znP)

        return xnP, znP, xnPQ, znPQ

    @staticmethod
    def _cubical_translate(xP, zP, xT, zT):
        r"""
        [Algorithm 4.103] CubicalTranslate, on coordinates.
        """
        x = xT * xP - zT * zP
        z = zT * xP - xT * zP

        if zT == 0:
            z = -z

        if xT == 0:
            x = -x

        return x, z

    def tate_pairing(self, Q, PQ, e):
        r"""
        [Algorithm 4.104] Tate.

        INPUT:

        - ``Q`` -- a point
        - ``PQ`` -- the point `P - Q`
        - ``e`` -- a positive integer with `[2^e]P = 0`

        OUTPUT:

        The reduced Tate pairing `t_{2^e}(P, Q)`.
        """

        p = self._parent._p

        xP, zP = self._x, self._z
        xP = xP / zP
        xQ, zQ = Q.xz()
        xQ = xQ / zQ
        xPQ, zPQ = PQ.xz()
        xPQ = xPQ / zPQ

        xnP, znP, xnPQ, znPQ = self._cubical_ladder(e-1, xPQ, 1, xP, 1, xQ)
        xnPQ, znPQ = self._cubical_translate(xnPQ, znPQ, xnP, znP)
        xnP, znP = self._cubical_translate(xnP, znP, xnP, znP)
        if znP != 0:
            raise ValueError("Tate failed: P does not have order 2^e")

        λ = znPQ / xnP
        return λ ** ((p**2 - 1) / (2**e))

    def has_order_2e(self, e):
        r"""
        Return whether the point has order exactly `2^e`.
        """
        P2 = self.ec_dbl_iter(e - 1)
        inf = P2.ec_dbl()

        _, Z2 = P2.xz()
        _, zero = inf.xz()

        return Z2 != 0 and zero == 0
