from sage.all import GF, is_prime, proof, log, valuation, matrix, ceil, Zmod

from ec.ec_curve import ECCurve
from ec.ec_basis import ECBasis

import params.p_324_3 as p324
import params.p_500_27 as p500
import params.p_664_17 as p664

proof.all(False)

class ECContext:

    def __init__(self, p, lb):
        r"""
        EllipticCurve context, derived from ``p`` and ``lb``.

        INPUT:

        - ``p`` -- the base prime
        - ``lb`` -- the security parameter `\lambda`
        """
        self.p = p
        self.lb = lb

        if not is_prime(p):
            raise ValueError(f"{p = } not prime")

        if p % 4 != 3:
            raise ValueError("p != 3 mod 4")

        self.Fp2, self._i = GF((self.p, 2), name='i', modulus=[1,0,1]).objgen()

        self.E0 = ECCurve(self.Fp2(0))
        self.FP_ENCODED_BYTES = ceil(log(self.p, 2) / 8)
        self.CHALLENGE_BYTES = lb // 8

        self.f = valuation(self.p+1, 2)
        self.e_rsp = 1 + ceil((lb + 2 * log(p, 2))/4)
        self.RESPONSE_BYTES = ceil((self.e_rsp + 2) / 8)

        # Load all precomputed points
        i = self._i
        if p == 3 * 2**324 - 1:
            self._P0 = self.E0([p324.xP0[0] + i*p324.xP0[1], 1])
            self._Q0 = self.E0([p324.xQ0[0] + i*p324.xQ0[1], 1])
            self._PQ0 = self.E0([p324.xPQ0[0] + i*p324.xPQ0[1], 1])
            self.det_basis_E0 = ECBasis(self._P0, self._Q0, self._PQ0)

            self.M_1 = matrix(Zmod(2**self.f), 2, [1, 0, 0, 1])
            self.M_i = matrix(Zmod(2**self.f), 2, [
                p324.aiP0, p324.aiQ0,
                p324.biP0, p324.biQ0
            ])

            self.M_ij = matrix(Zmod(2**self.f), 2, [
                p324.aijP0, p324.aijQ0,
                p324.bijP0, p324.bijQ0
            ])

            self.M_1k = matrix(Zmod(2**self.f), 2, [
                p324.a1kP0, p324.a1kQ0,
                p324.b1kP0, p324.b1kQ0
            ])

            self.matrices = [self.M_1, self.M_i, self.M_ij, self.M_1k]

        elif p == 27 * 2**500 - 1:
            self._P0 = self.E0([p500.xP0[0] + i*p500.xP0[1], 1])
            self._Q0 = self.E0([p500.xQ0[0] + i*p500.xQ0[1], 1])
            self._PQ0 = self.E0([p500.xPQ0[0] + i*p500.xPQ0[1], 1])
            self.det_basis_E0 = ECBasis(self._P0, self._Q0, self._PQ0)

            self.M_1 = matrix(Zmod(2**self.f), 2, [1, 0, 0, 1])
            self.M_i = matrix(Zmod(2**self.f), 2, [
                p500.aiP0, p500.aiQ0,
                p500.biP0, p500.biQ0
            ])

            self.M_ij = matrix(Zmod(2**self.f), 2, [
                p500.aijP0, p500.aijQ0,
                p500.bijP0, p500.bijQ0
            ])

            self.M_1k = matrix(Zmod(2**self.f), 2, [
                p500.a1kP0, p500.a1kQ0,
                p500.b1kP0, p500.b1kQ0
            ])

            self.matrices = [self.M_1, self.M_i, self.M_ij, self.M_1k]

        elif p == 17 * 2**664 - 1:
            self._P0 = self.E0([p664.xP0[0] + i*p664.xP0[1], 1])
            self._Q0 = self.E0([p664.xQ0[0] + i*p664.xQ0[1], 1])
            self._PQ0 = self.E0([p664.xPQ0[0] + i*p664.xPQ0[1], 1])
            self.det_basis_E0 = ECBasis(self._P0, self._Q0, self._PQ0)

            self.M_1 = matrix(Zmod(2**self.f), 2, [1, 0, 0, 1])
            self.M_i = matrix(Zmod(2**self.f), 2, [
                p664.aiP0, p664.aiQ0,
                p664.biP0, p664.biQ0
            ])

            self.M_ij = matrix(Zmod(2**self.f), 2, [
                p664.aijP0, p664.aijQ0,
                p664.bijP0, p664.bijQ0
            ])

            self.M_1k = matrix(Zmod(2**self.f), 2, [
                p664.a1kP0, p664.a1kQ0,
                p664.b1kP0, p664.b1kQ0
            ])

            self.matrices = [self.M_1, self.M_i, self.M_ij, self.M_1k]

        else:
            raise NotImplementedError('prime not implemented')
