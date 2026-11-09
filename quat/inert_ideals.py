from sage.all import (
        ZZ, Zmod, crt, matrix, QQ, pari, QuadraticField, isqrt, vector, floor
)

import quat.quaternions as qt
import crypto.bint as bint
import crypto.utils as utils

def ideal_to_Nxy(I):
    r"""
    Return the inert form of ``I``. Inverse of :func:`xy_to_ideal`.

    INPUT:

    - ``I`` -- an inert `O_0`-ideal

    OUTPUT:

    A tuple ``(N, x, y)`` of integers with ``N`` the norm of ``I``.
    """
    # Ideal HNF
    C = matrix(QQ, [qt.coordinates_O0(b) for b in I.basis()])
    H = pari.mathnf(C.change_ring(ZZ).transpose()).sage().transpose()

    N, x, y = H[0, 0], H[2, 0], H[2, 1]

    p = I.quaternion_algebra().ramified_primes()[0]
    if (H != matrix(ZZ, [[N, 0, 0, 0], [0, N, 0, 0],
                         [x, y, 1, 0], [(-y-1) % N, x, 0, 1]])
            or (4 * x**2 + (2 * y + 1) ** 2 + p) % (4 * N) != 0):
        raise ValueError("the ideal is not inert")

    return N, x, y

def xy_to_ideal(N, x, y, quat_ctx):
    r"""
    [Algorithm 4.66] InertLatticeBasis.

    Return the ideal with inert form ``(N, (x, y))``.

    INPUT:

    - ``N`` -- the ideal norm
    - ``x``, ``y`` -- integers in `[0, N)`
    - ``quat_ctx`` -- the quaternion context

    OUTPUT:

    The left `O_0`-ideal with basis
    `(N, Ni, x + yi + (i+j)/2, -y + xi + (k-1)/2)`.
    """
    i, j, k = quat_ctx.B.gens()
    id_basis = [N, N * i, x + y * i + (i + j) / 2, -y + x * i + (k - 1) / 2]
    return quat_ctx.O0.left_ideal(id_basis, is_basis=True)

def prime_norm_inert_ideal_sampling(N, quat_ctx, prng_key, prng_def):
    r"""
    [Algorithm 4.71] PrimeNormIdealInertSampling.

    Sample a random inert `O_0`-ideal of prime norm ``N``.

    INPUT:

    - ``N`` -- an odd prime
    - ``quat_ctx`` -- the quaternion context
    - ``prng_key`` -- the pseudorandom generator for sampling
    - ``prng_def`` -- the default-domain pseudorandom generator

    OUTPUT:

    A tuple ``(x, y)`` such that ``(N, x, y)`` is the inert form of the ideal.
    """
    while True:
        # Here the subtraction by 1 increases the bitbound on N by one; it's
        # therefore important to cast before calling the prng
        x = int(prng_key.sample_from_interval(0, bint.BoundedInt(N)-1))
        t = -4*x**2 - quat_ctx.p
        if Zmod(4*N)(t).is_square():
            r = utils.modular_sqrt(t, N, prng_def)
            r4 = t % 2
            Y = crt(r4, r, 4, N)
            b = prng_key.sample_from_interval(0, 1)
            Y = ((-1)**b * Y) % (2*N)
            y = (Y - 1)/2
            return x, y

def inert_ideal_reduction(I, p):
    r"""
    Replaces [Algorithm 4.72] InertIdealReduction.

    Return a canonical reduced basis of an inert ideal.

    INPUT:

    - ``I`` -- an inert ideal ``(N, x, y)``
    - ``p`` -- the base prime

    OUTPUT:

    A tuple ``(a, c)`` of coordinate lists in the basis `(1, i, j, k)`, such
    that `(a, ia, c, ic)` is a reduced basis of ``I``.
    """
    K, ii = QuadraticField(-1, names='ii').objgen()

    N, x, y = I
    X, Y, qh = 2*x, 2*y + 1, isqrt(p)
    v = vector(K, [2*N, 0])
    w = vector(K, [X + Y*ii, qh])

    def H(a, b):
        # Hermitian prod
        return a[0].conjugate()*b[0] + a[1].conjugate()*b[1]

    def rnd(z):
        # Nearest gaussian integer
        return floor(z[0] + 1/2) + floor(z[1] + 1/2)*ii

    # Euclidean algorithm
    while True:
        if H(w, w) < H(v, v):
            v, w = w, v
        r = rnd(H(v, w) / H(v, v))
        if r == 0:
            break
        w = w - r*v

    def orient_unit(z):
        # The unit u such that Re(u*z) >= 0 and Im(u*z) >= 0
        u = K(1)
        if (z[0] < 0) != (z[1] < 0):
            u, z = ii, ii*z
        if z[0] < 0 or z[1] < 0:
            u, z = -u, -z
        return u

    # Put v and <v, w> in the first quadrant
    v = orient_unit(v[0]) * v
    h = v[0].conjugate()*w[0] + p * (v[1]/qh).conjugate() * (w[1]/qh)
    w = orient_unit(h) * w

    # Undo the sqrt(p) scaling
    def zi_to_quat(z):
        z1, z2 = z[0], z[1] / qh
        return [z1[0] / 2, z1[1] / 2, z2[0] / 2, z2[1] / 2]

    a = zi_to_quat(v)
    c = zi_to_quat(w)
    return a, c
