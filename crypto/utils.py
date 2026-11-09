from sage.all import ZZ, Zmod, kronecker

import crypto.bint as bint
import quat.quaternions as qt

def modular_sqrt(x, p, prng):
    r"""
    [Algorithm 4.20] ModularSQRT.

    INPUT:

    - ``x`` -- an integer
    - ``p`` -- an odd prime, as an integer or :class:`BoundedInt`
    - ``prng`` -- the pseudorandom generator, used for Tonelli-Shanks when
      `p \equiv 1 \bmod 4`

    OUTPUT:

    An integer ``r`` in `[0, p)` with `r^2 \equiv x \bmod p`, or ``None`` if
    ``x`` is not a square.
    """
    if not isinstance(p, bint.BoundedInt):
        p = bint.BoundedInt(p)
    pp = ZZ(p.nn)
    F = Zmod(pp)
    a = F(ZZ(int(x)))
    if a == 0:
        return ZZ(0)

    if pp % 4 == 3:
        # p = 3 mod 4
        r = a**((pp + 1) // 4)
        return ZZ(r) if r**2 == a else None
    if pp % 4 != 1:
        return None

    # p = 1 mod 4: Tonelli-Shanks
    pm1 = p - 1
    s = ZZ(pm1.nn).valuation(2)
    q = ZZ(pm1.nn) >> s

    # Quadratic non-residue z sampled in [1, p-1]
    while True:
        z = ZZ(prng.sample_from_interval(1, pm1))
        if kronecker(z, pp) == -1:
            break

    c = F(z)**q
    t = a**q
    R = a**((q + 1) // 2)

    M = s
    while t != 1:
        tmp = t
        i = 0
        while tmp != 1:
            tmp = tmp**2
            i += 1
            if i >= M:
                return None

        b = c**(2**(M - i - 1))
        R = R * b
        c = b**2
        t = t * c
        M = i

    return ZZ(R)

def fp2_sqrt(x):
    r"""
    [Algorithm 4.35] Fp2SQRT.

    INPUT:

    - ``x`` -- a square in `\mathbb{F}_{p^2}`

    OUTPUT:

    An element ``r`` with `r^2 = x`.
    """
    Fp2 = x.parent()
    p = Fp2.characteristic()

    x0, x1 = x
    delta = (x0**2 + x1**2) ** ((p+1)//4)

    if x1 == 0:
        delta = x0

    r0 = x0 + delta
    t = 2*r0
    r1 = t**((p-3)//4)
    r0 = r0 * r1
    r1 = x1 * r1

    if (2 * r0)**2 == t:
        return Fp2([r0, r1])
    return Fp2([r1, -r0])

def normalized_dlog(f, g_inv, e):
    r"""
    [Algorithm 4.37] NormalizedDlog.

    INPUT:

    - ``f`` -- an element of `\langle g \rangle`
    - ``g_inv`` -- the inverse of an element `g` of order `2^e`
    - ``e`` -- a non-negative integer

    OUTPUT:

    The integer ``a`` in `[0, 2^e)` with `f = g^a`.
    """
    if e == 0:
      return 0
    if e == 1:
      if f == 1:
          return 0
      return 1

    r = e // 2
    l = e - r

    a0 = normalized_dlog(f**(2**l), g_inv**(2**l), r)
    f = f * g_inv**a0
    g_inv = g_inv**(2**r)

    a1 = normalized_dlog(f, g_inv, l)
    return a0 + 2**r * a1

def fp2_compare(a, b):
    r"""
    Compare two elements of `\mathbb{F}_{p^2}` [Definition 4.2.1].

    INPUT:

    - ``a``, ``b`` -- elements of `\mathbb{F}_{p^2}`

    OUTPUT:

    ``True`` if `a < b`, ``False`` otherwise.
    """
    return (a[1] < b[1]) or (a[1] == b[1] and a[0] < b[0])

def batched_inversion(*input):
    r"""
    [Algorithm 4.36] Fp2BatchInverse.

    INPUT:

    - ``input`` -- invertible elements `x_1, \ldots, x_n`

    OUTPUT:

    The list `[x_1^{-1}, \ldots, x_n^{-1}]`.
    """
    multiples = [input[0]]
    for ai in input[1:]:
        multiples.append(ai * multiples[-1])

    last_multiple = multiples[-1]
    inverses_multiples = [1 / last_multiple]

    for ai in reversed(input[1:]):
        inverses_multiples.append(ai * inverses_multiples[-1])

    inverses_multiples = inverses_multiples[::-1]
    inverses = [inverses_multiples[0]]

    k = len(input)
    for i in range(1, k):
        inverses.append(inverses_multiples[i] * multiples[i - 1])

    return inverses

def set_change_of_basis_matrix(E1, E2, B1, B2, e):
    r"""
    [Algorithm 3.9] SetChangeOfBasisMatrix.

    INPUT:

    - ``E1``, ``E2`` -- curves of type :class:`ECCurve`
    - ``B1``, ``B2`` -- bases of `E_1[2^e]` and `E_2[2^e]`
    - ``e`` -- an integer with `e \le f`

    OUTPUT:

    A tuple ``(M2, hint_1, hint_2)`` with ``M2`` a `2 \times 2` matrix over
    `\ZZ/2^e\ZZ` and the hints of the canonical bases of ``E1`` and ``E2``.
    """
    f = (E1._p + 1).valuation(2)
    B1_1, hint_1 = E1.torsion_basis_to_hint(f, e)
    B2_1, hint_2 = E2.torsion_basis_to_hint(f, e)

    M1 = B1_1.change_of_basis(B1, e)

    B2 = B2.act_on_basis(M1**-1)
    M2 = B2_1.change_of_basis(B2, e)

    return M2, hint_1, hint_2

def normalize_matrix(M, e):
    r"""
    [Algorithm 3.10] NormalizeMatrix.

    INPUT:

    - ``M`` -- a tuple of four integers in `[0, 2^e)`
    - ``e`` -- a positive integer

    OUTPUT:

    The lexicographically smallest among ``M`` and `-M \bmod 2^e`.
    """
    for i in range(4):
        if 0 < ZZ(M[i]) < 2**(e-1):
            return M
        elif ZZ(M[i]) > 2**(e-1):
            return tuple((-x % 2**e) for x in M)
    return M

def endomorphism_to_matrix(ec_ctx, beta):
    r"""
    [Algorithm 4.134] EndomorphismToMatrix.

    INPUT:

    - ``ec_ctx`` -- the elliptic curve context
    - ``beta`` -- an element of `O_0`

    OUTPUT:

    The matrix over `\ZZ/2^f\ZZ` of the action of ``beta`` on the basis
    `(P_0, Q_0)` of `E_0[2^f]`.
    """
    x, y, z, t = map(ZZ, qt.coordinates_O0(beta))
    mat_1, mat_i, mat_ip2, mat_1ip2 = ec_ctx.matrices

    return x*mat_1 + y*mat_i + z*mat_ip2 + t*mat_1ip2
