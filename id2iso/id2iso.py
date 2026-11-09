from sage.all import (
    gcd, inverse_mod, is_prime, log, matrix, randint, sqrt, Zmod, ZZ
)

import quat.quaternions as qt
import quat.inert_ideals as qt_inert
import id2iso.dim_2_lattice as d2l

from hd.couple_point import CoupleECPoint
from hd.isogeny_chain import isogeny_two_two_chain

from ec.ec_basis import ECBasis
import crypto.utils as utils


def ideal_to_isogeny(J, prng, quat_ctx, ec_ctx):
    r"""
    [Algorithm 4.150] IdealToIsogeny.

    Return the codomain of `\varphi_J` and the image of the torsion basis of
    `E_0`.

    INPUT:

    - ``J`` -- a left `O_0`-ideal
    - ``prng`` -- the pseudorandom generator
    - ``quat_ctx`` -- the quaternion context
    - ``ec_ctx`` -- the elliptic curve context

    OUTPUT:

    A tuple ``(E_I, B_I)`` with ``E_I`` the codomain of `\varphi_J` and
    ``B_I`` the basis `(\varphi_J(P_0), \varphi_J(Q_0), \varphi_J(P_0 - Q_0))`.
    """
    f = ec_ctx.f
    B0 = ec_ctx.det_basis_E0
    P0, Q0, _ = B0

    beta1, d_1, theta = ideal_to_isogeny_norm_equation(quat_ctx, J, f - 2, prng)

    d_2 = 2**(f - 2) - ZZ(d_1)
    modulus = 2**f
    nu_1 = inverse_mod(d_1, modulus)
    nu_2 = inverse_mod(d_2, modulus)
    c_2 = nu_1 * (1 - nu_2 * 2**(f - 2)) % modulus

    M_theta = utils.endomorphism_to_matrix(ec_ctx, theta)
    M_theta *= matrix(Zmod(2**f), 2, [nu_1, 0, 0, c_2])
    T12, T22, _ = B0.act_on_basis(M_theta)

    T1 = CoupleECPoint(P0, T12)
    T2 = CoupleECPoint(Q0, T22)

    M_beta1 = utils.endomorphism_to_matrix(ec_ctx, beta1)
    P1, Q1, R1 = B0.act_on_basis(nu_1 * M_beta1)

    zero = ec_ctx.E0.zero()
    U1 = CoupleECPoint(P1, zero)
    V1 = CoupleECPoint(Q1, zero)
    W1 = CoupleECPoint(R1, zero)

    E_I, _, image_points = isogeny_two_two_chain(
        ec_ctx.E0, ec_ctx.E0, T1, T2, f-2, [U1, V1, W1], False)
    U_I, V_I, W_I = image_points
    P_I, _ = U_I.points()
    Q_I, _ = V_I.points()
    R_I, _ = W_I.points()

    return E_I, ECBasis(P_I, Q_I, R_I)

def gauss_gcd(quat_ctx, alpha, beta):
    r"""
    [Algorithm 4.48] GaussGCD.

    INPUT:

    - ``quat_ctx`` -- the quaternion context
    - ``alpha``, ``beta`` -- elements of `\ZZ[i]`

    OUTPUT:

    A gcd of ``alpha`` and ``beta`` in `\ZZ[i]`.
    """
    kappa, rho = alpha, beta

    if kappa.reduced_norm() < rho.reduced_norm():
        kappa, rho = rho, kappa

    if rho.reduced_norm() == 1:
        return rho

    if rho == 0:
        return kappa

    a, b = list(kappa/rho)[:2]
    gamma = a.round("away") + quat_ctx.i*b.round("away")
    delta = kappa - gamma*rho

    return gauss_gcd(quat_ctx, delta, rho)

def loop_one(quat_ctx, N, beta_I, a, b, v0, prng):
    r"""
    [Algorithm 4.146] LoopOne.

    INPUT:

    - ``quat_ctx`` -- the quaternion context
    - ``N`` -- a positive integer
    - ``beta_I`` -- an element of `O_0`
    - ``a``, ``b`` -- bounds with `0 \le a < b`
    - ``v0`` -- the smallest vector in the kernel lattice of ``beta_I``
    - ``prng`` -- the pseudorandom generator

    OUTPUT:

    A tuple ``(alpha0, Lambda, v1)`` with ``alpha0`` in
    `O_0 \langle \beta_I, N \rangle`, ``Lambda`` the reduced kernel lattice
    of ``alpha0`` and ``v1`` a 2D vector mapped to `1 \bmod 4N`.
    """
    s0, t0 = v0
    ell = gcd(s0, t0)
    omega = s0 + quat_ctx.i*t0
    N_four = 4*N
    trace_pairing = lambda alpha, v: \
        ((v[0] - quat_ctx.i*v[1]) * alpha).reduced_trace() % N_four

    if trace_pairing(beta_I, v0) != 0:
        raise ValueError("v0 is not in the kernel lattice of beta_I")

    while True:
        if b < 5:
            c3 = 1
            c4 = 0
        else:
            c3 = randint(sqrt(a/2).n().ceil(), sqrt(b).n().floor())
            c4 = randint(sqrt(max(0, a-c3**2)).n().ceil(), sqrt(b - c3**2).n().floor())
        if c3 % 2 == c4 % 2 or gcd(c3, c4) != 1:
            continue
        if sqrt(b) - sqrt(a) >= log(N)**2:
            if gcd(c3**2 + c4**2, N) != 1:
                continue
        else:
            if gcd(c3**2 + c4**2, ell) != 1:
                continue
        z_c = c3 + quat_ctx.i*c4
        gamma = gauss_gcd(quat_ctx, z_c.conjugate(), omega)

        while gamma.reduced_norm() != 1:
            z_c = z_c*gamma**2/gamma.reduced_norm()
            gamma = gauss_gcd(quat_ctx, z_c.conjugate(), omega)

        alpha0 = z_c * beta_I
        w = list(z_c * (s0 + quat_ctx.i*t0))[:2]

        a, b, c, d, e = qt.coordinates_B(alpha0)

        a = a % (N_four)
        b = b % (N_four)
        alpha0 = (a + quat_ctx.i*b + quat_ctx.j*c + quat_ctx.k*d)/2

        Lambda, v1 = d2l.compute_2d_lattice(quat_ctx, N, alpha0)

        Lambda = d2l.reduce_given_minima(Lambda, w)

        return alpha0, Lambda, v1


def loop_two(quat_ctx, e, N, alpha0, Lambda, v):
    r"""
    [Algorithm 4.148] LoopTwo.

    INPUT:

    - ``quat_ctx`` -- the quaternion context
    - ``e`` -- a positive integer
    - ``N`` -- a positive integer
    - ``alpha0`` -- an element of `O_0` with `\mathrm{nrd}(\alpha_0) \equiv
      0 \bmod N`
    - ``Lambda`` -- the kernel lattice of ``alpha0``
    - ``v`` -- a 2D vector mapped to `1 \bmod 4N`

    OUTPUT:

    A tuple ``(beta1, beta2)`` of elements of `O_0 \langle \alpha_0, N
    \rangle` with `\mathrm{nrd}(\beta_1) + \mathrm{nrd}(\beta_2) = N 2^e`,
    or ``None`` if no solution is found.
    """
    r0 = alpha0.reduced_norm()//N
    delta = 1
    if r0 % 2 == 0:
        delta = 3
    lam = -1

    lam_bound = sqrt(2**(e-2)/r0).n()
    while lam < lam_bound:
        lam = lam + 2
        if gcd(2*N, lam) != 1:
            continue

        d = 2**e - 2*lam**2*r0 - delta*N
        s, t = d2l.compute_small_solutions(N, Lambda, lam,  v, d)

        if s % 2 == t % 2:
            continue

        a = alpha0.reduced_trace()
        b = -(quat_ctx.i*alpha0).reduced_trace()
        k = ZZ((2**e - 2*lam**2*r0 - lam*a*s - lam*b*t)/N)
        z = 2*k - s**2 - t**2

        if z < 0:
            continue

        if not is_prime(z):
            continue
        sol = qt.cornacchia(z)
        if not sol:
            continue

        sc, tc = sol
        if sc % 2 != s % 2:
            sc, tc = tc, sc

        a1 = (sc + s)/2
        b1 = (tc + t)/2
        a2 = s - a1
        b2 = t - b1

        beta1 = N*(a1 + quat_ctx.i*b1) + lam*alpha0
        beta2 = N*(a2 + quat_ctx.i*b2) + lam*alpha0

        return beta1, beta2
    return None


def ideal_to_isogeny_norm_equation(quat_ctx, J, e, prng):
    r"""
    [Algorithm 4.149] IdealToIsogenyNormEquation.

    INPUT:

    - ``quat_ctx`` -- the quaternion context
    - ``J`` -- a left `O_0`-ideal
    - ``e`` -- a positive integer
    - ``prng`` -- the pseudorandom generator

    OUTPUT:

    A tuple ``(beta1, d1, theta)`` with ``beta1`` in ``J``,
    `d_1 = \mathrm{nrd}(\beta_1)/\mathrm{nrd}(J)` and
    `\theta = \beta_2 \bar\beta_1 / \mathrm{nrd}(J)`, where
    `\mathrm{nrd}(\beta_1) + \mathrm{nrd}(\beta_2) = 2^e \mathrm{nrd}(J)`.
    """
    gamma = J.reduced_basis()[0]

    I = J * (gamma.conjugate() / J.norm())
    N, x, y = qt_inert.ideal_to_Nxy(I)

    B_1 = (N/2).n()
    B_M = ((5*N)/2).n()
    B_alpha = (N*(2**(e+3) - 16*N)/quat_ctx.p).n()

    beta_I = x + quat_ctx.i*y + (quat_ctx.i + quat_ctx.j)/2

    Lambda, _ = d2l.compute_2d_lattice(quat_ctx, N, beta_I)
    Lambda = d2l.lagrange_reduction(Lambda)
    v0, _ = d2l.lat_to_basis(Lambda)
    lambda1_m = d2l.vec_norm(v0)

    a = (B_1/lambda1_m).n()
    b = min(B_alpha/2, B_M/lambda1_m).n()

    if sqrt(b) - sqrt(a) <= log(N)**2:
        beta_I, a, b, v0 = d2l.divide_minima(quat_ctx, N, beta_I, a, b, v0)

    found = False
    while not found:
        alpha0, Lambda, v1 = loop_one(quat_ctx, N, beta_I, a, b, v0, prng)
        sol = loop_two(quat_ctx, e, N, alpha0, Lambda, v1)
        if sol:
            found = True

    beta1, beta2 = sol
    if beta1.reduced_norm() + beta2.reduced_norm() != N*2**e:
        raise RuntimeError("IdealToIsogenyNormEquation failed: wrong norms")
    d1 = ZZ(beta1.reduced_norm()/N)
    theta = (beta2*beta1.conjugate())/N
    beta1 = (beta1*gamma)/N

    return beta1, d1, theta
