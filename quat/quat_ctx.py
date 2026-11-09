from sage.all import (
    valuation,
    QuaternionAlgebra,
    is_prime,
    proof,
    next_prime,
    ZZ,
    floor,
    ceil,
    log,
)

proof.all(False)

def compute_EIBox(p, lb):
    r"""
    [Algorithm B.2] EIBox.

    Return the least integer `b` such that `(1 - 1/(2\log p))^{(b+1)^4} < 2^{-\lambda}`.

    INPUT:

    - ``p`` -- the base prime
    - ``lb`` -- the security parameter `\lambda`

    OUTPUT:

    An integer ``b``.
    """
    b = ceil(log(p, 2))
    b = log(1- 1/(2*b), 2)
    b = -lb/b
    b = b**(1/4) - 1
    b = ceil(b)
    return b

def compute_RICofactor(p, lb):
    r"""
    [Algorithm B.3] RICofactor.

    Return the cofactor used in the norm of RepresentInteger.

    INPUT:

    - ``p`` -- the base prime
    - ``lb`` -- the security parameter `\lambda`

    OUTPUT:

    A prime ``q``.
    """
    L = lambda n, b: round(float(log(1 - (1/(2*b)), 2) * n))

    e_rsp = 1 + ceil((lb + 2 * log(p, 2))/4)
    m = ceil(log(p, 2)) - e_rsp + 2
    while L(floor(2**(e_rsp+m) / p), e_rsp + m) > -lb:
        m += 1
    q = next_prime(2**m)
    return q

class QuatContext:
    def __init__(self, p, lb):
        r"""
        Quaternion context, derived from ``p`` and ``lb``.

        INPUT:

        - ``p`` -- the base prime, `p \equiv 3 \bmod 4`
        - ``lb`` -- the security parameter `\lambda`
        """
        self.p = ZZ(p)

        if not is_prime(p):
            raise ValueError(f"{p = } not prime")

        if p % 4 != 3:
            raise ValueError("p != 3 mod 4")

        self.B = QuaternionAlgebra(-1, -p)
        i, j, k = self.B.gens()
        self.i = i
        self.j = j
        self.k = k
        self.f = valuation(self.p + 1, 2)
        self.O0_basis = [1, i, (i + j) / 2, (1 + k) / 2]
        self.O0 = self.B.quaternion_order(self.O0_basis)

        self.m = next_prime(2 ** (self.p.nbits()))
        self.Dmix = next_prime(self.p * 2 ** (2 * lb))
        self.EIBox = compute_EIBox(p, lb)
        self.EIBoundIteration = (self.EIBox + 1)**4
        self.RICofactor = compute_RICofactor(p, lb)
        self.BoundIterationRBE = ceil(11.636 * lb)
