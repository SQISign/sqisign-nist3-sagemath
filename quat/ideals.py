from sage.all import (
        matrix, gcd, QQ, ZZ, is_prime, isqrt, vector, Zmod, QuadraticForm,
        GaussianIntegers
)

import quat.quaternions as qt
import quat.inert_ideals as qt_inert
import crypto.bint as bint

def gen_to_ideal(N, gamma, quat_ctx):
    r"""
    [Algorithm 4.76] GenToIdealOdd.

    Return the ideal `O_0 \langle \gamma, N \rangle` in generic form.

    INPUT:

    - ``N`` -- a positive odd integer coprime to `p`
    - ``gamma`` -- an element of `O_0` with `\mathrm{nrd}(\gamma) \equiv 0 \bmod N`
    - ``quat_ctx`` -- the quaternion context

    OUTPUT:

    A tuple ``(I, theta)`` with ``I`` an ideal and ``theta = u + v*i`` such
    that `O_0 \langle \gamma, N \rangle = I \theta`.
    """
    ii = GaussianIntegers().gen(1)
    gamma_N = qt.mod_O0(gamma, N)
    gamma_N = qt.from_O0_coords(gamma_N, quat_ctx)
    a, b, c, d, e = qt.coordinates_B(gamma_N)
    uv = gcd(a + b*ii, c - d*ii)
    uv = gcd(uv, N)

    if uv in [ii, -ii]:
        uv = (1, 0)
    u, v = uv

    # Return the ideal
    i = quat_ctx.B.gens()[0]
    O0 = quat_ctx.O0

    N2 = u**2 + v**2
    N1 = N / N2
    gamma *= (u - i*v)/N2
    return O0*N1 + O0*gamma, u + i*v

def random_ideal_given_norm(N, alpha, prng, quat_ctx):
    r"""
    [Algorithm 4.82] RandomIdealGivenNorm.

    Sample a random `O_0`-ideal of norm ``N``.

    INPUT:

    - ``N`` -- a positive odd integer coprime to `p`
    - ``alpha`` -- an element of `O_0` of norm coprime to ``N``
    - ``prng`` -- the pseudorandom generator
    - ``quat_ctx`` -- the quaternion context

    OUTPUT:

    A tuple ``(I, theta)`` as in :func:`gen_to_ideal`.
    """
    i, j, k = quat_ctx.B.gens()
    gamma = qt.represent_integer(N * quat_ctx.RICofactor, prng, quat_ctx)

    while True:
        x, y, z, w = [prng.sample_from_interval(1, N) for _ in range(4)]
        beta = x + y*i + z*j + w*k
        M = beta.reduced_norm() % N
        if gcd(M, N) == 1:
            break

    gamma1 = gamma * beta * alpha
    I, gamma = gen_to_ideal(N, gamma1, quat_ctx)
    return I, gamma

def sample_random_from_basis(I_basis, b, prng):
    r"""
    [Algorithm 4.78] SampleRandomFromBasis.

    Sample a primitive element of the lattice spanned by ``I_basis``, with
    coefficients in `[-b, b]`.

    INPUT:

    - ``I_basis`` -- a list of four quaternion elements
    - ``b`` -- a positive integer
    - ``prng`` -- the pseudorandom generator

    OUTPUT:

    A quaternion element ``beta``.
    """
    # This is matching the C code since the spec differ
    mm = bint.BoundedInt(b, 32)
    mm = mm + mm
    coeffs = [prng.sample_from_interval(0, mm) - b for _ in range(4)]

    # Make primitive
    g = gcd(coeffs)
    coeffs = [c // g for c in coeffs]

    beta = sum(coeffs[i] * I_basis[i] for i in range(4))
    return beta

def equivalent_coprime_ideal(I, M, prng, quat_ctx):
    r"""
    [Algorithm 4.79] EquivalentCoprimeIdeal.

    Return an ideal equivalent to ``I`` of norm coprime to ``M`` if
    `M \neq 0`, and of prime norm otherwise.

    INPUT:

    - ``I`` -- an inert ideal ``(N, x, y)``
    - ``M`` -- an integer
    - ``prng`` -- the pseudorandom generator
    - ``quat_ctx`` -- the quaternion context

    OUTPUT:

    A tuple ``(J, beta)`` with ``J`` the equivalent ideal and ``beta`` the
    element of ``I`` such that `J = I \bar\beta / \mathrm{nrd}(I)`.
    """
    B = quat_ctx.B
    i, j, k = B.gens()

    a, c = qt_inert.inert_ideal_reduction(I, quat_ctx.p)
    I_basis = [B(a), i*B(a), B(c), i*B(c)]
    N, x, y = I
    I = qt_inert.xy_to_ideal(N, x, y, quat_ctx)

    for _ in range(quat_ctx.EIBoundIteration):
        beta = sample_random_from_basis(I_basis, quat_ctx.EIBox, prng)
        N1 = ZZ(beta.reduced_norm() / I.norm())
        if M != 0:
            # Sample ideal of norm coprime with M for signing
            if gcd(M, N1) == 1:
                J = I * beta.conjugate() * (1 / I.norm())
                return J, beta
        else:
            # M = 0 -> sample prime norm
            if is_prime(N1):
                J = I * beta.conjugate() * (1 / I.norm())
                return J, beta

    raise RuntimeError("equivalent_coprime_ideal failed")

def uniform_O0_class_sampling(M, prng_key, prng_def, quat_ctx):
    r"""
    [Algorithm 4.81] UniformO0ClassSampling.

    Sample a random `O_0`-ideal class.

    INPUT:

    - ``M`` -- an integer
    - ``prng_key`` -- the pseudorandom generator for sampling
    - ``prng_def`` -- the default-domain pseudorandom generator
    - ``quat_ctx`` -- the quaternion context

    OUTPUT:

    An ideal ``I`` of norm coprime to ``M`` if `M \neq 0`, and prime
    otherwise.
    """
    x, y = qt_inert.prime_norm_inert_ideal_sampling(quat_ctx.Dmix, quat_ctx,
                                                prng_key, prng_def)
    I = (quat_ctx.Dmix, x, y)
    I, _ = equivalent_coprime_ideal(I, M, prng_key, quat_ctx)

    return I

def kernel_to_ideal(e, c1, c2, ec_ctx, quat_ctx):
    r"""
    [Algorithm 4.151] KernelDecomposedToIdeal.

    Return the ideal of the isogeny with kernel `\langle [c_1]P_0 + [c_2]Q_0
    \rangle` on `E_0`.

    INPUT:

    - ``e`` -- a positive integer
    - ``c1``, ``c2`` -- integers defining a point of order `2^e`
    - ``ec_ctx`` -- the elliptic curve context
    - ``quat_ctx`` -- the quaternion context

    OUTPUT:

    A tuple ``(I, gamma)`` such that `I \gamma` is the ideal of norm `2^e`.
    """
    M_j = 2*ec_ctx.M_ij - ec_ctx.M_i
    M_alpha = M_j + ec_ctx.M_1k
    M_alpha = M_alpha.change_ring(Zmod(2**e))
    vc = vector(Zmod(2**e), [c1, c2])
    d1, d2 = M_alpha * vc
    M = matrix(Zmod(2**e), [[c1, d1], [c2, d2]])
    a, b = (M**-1 * ec_ctx.M_i.change_ring(Zmod(2**e))) * vc

    B = quat_ctx.B
    i, j, k = B.gens()
    alpha = ZZ(a) + ZZ(b)*(j + (1+k)/2) - i

    alpha1 = alpha * (1 + i)
    t, x, y, z = qt.coordinates_O0(alpha1)
    gamma = 1

    if all([cff % 2 == 0 for cff in [t, x, y, z]]):
        alpha = alpha1 / 2
        e -= 1
        gamma = 1 - i

    return quat_ctx.O0 * alpha + quat_ctx.O0 * 2**e, gamma


def dual_reduction(I):
    r"""
    Replaces [Algorithm 4.75] MulIdealDualReduction.

    Minkowski-reduce the dual lattice of ``I``.

    INPUT:

    - ``I`` -- an ideal

    OUTPUT:

    A tuple ``(lam, T)`` with ``lam`` the squared successive minima of the
    dual of ``I`` scaled by `\mathrm{nrd}(I)`, and ``T`` a `4 \times 4`
    integer matrix as in [Algorithm 4.75].
    """

    def trace_pairing(x, y):
        return (x * y.conjugate()).reduced_trace() / 2

    def gram(basis):
        return matrix(QQ, 4, [[trace_pairing(x, y) for y in basis] for x in basis])

    def dual_basis(basis):
        return [sum(c * b for c, b in zip(row, basis))
                for row in gram(basis).inverse().rows()]

    def minkowski_basis(basis):
        G = gram(basis)
        G = (G.denominator() * G).change_ring(ZZ)
        U = G.LLL_gram()

        while True:
            U = matrix(sorted(U.columns(), key=lambda u: u * G * u)).transpose()
            # minkowski_reduction does not sort the minima
            M = QuadraticForm(ZZ, 2 * U.transpose() * G * U).minkowski_reduction()[1]
            if M == 1:
                break
            U = U * M

        v = [sum(c * b for c, b in zip(u, basis)) for u in U.columns()]

        # Canonical choice
        for i in range(1, 4):
            if trace_pairing(v[0], v[i]) < 0:
                v[i] = -v[i]
        return v

    N = I.norm()

    dual_bs = dual_basis(I.basis())
    dual_min = minkowski_basis(dual_bs)
    bprime = dual_basis(dual_min)[::-1]
    lam = [N * di.reduced_norm() for di in dual_min[::-1]]

    C = matrix(ZZ, 4, 4, lambda r, c: trace_pairing(bprime[c], dual_bs[r]))
    return lam, C[::-1, :]

def random_bounded_elt(I, B, prng, quat_ctx):
    r"""
    [Algorithm 4.83] RandomBoundedElement.

    Sample an element of ``I`` of small odd norm.

    INPUT:

    - ``I`` -- an ideal
    - ``B`` -- a positive integer bound
    - ``prng`` -- the pseudorandom generator
    - ``quat_ctx`` -- the quaternion context

    OUTPUT:

    A tuple ``(M, alpha)`` with ``alpha`` in ``I`` and `M =
    \mathrm{nrd}(\alpha)/\mathrm{nrd}(I)` odd and at most ``B``.
    """
    bs = I.basis()
    lambdas, T = dual_reduction(I)
    N = I.norm()

    betas = []
    for i in range(4):
        # Here we avoid the scaling since we did not scale before
        cff_i = B * lambdas[i]
        beta_i = bint.BoundedInt(isqrt(cff_i))
        betas.append(beta_i)

    for _ in range(quat_ctx.BoundIterationRBE):
        xi = [prng.sample_from_interval(0, betas[i] + betas[i]) - betas[i].nn
              for i in range(4)]
        Xi = vector(ZZ, xi)

        y = (T * Xi)[::-1]
        alpha = sum(yi * bi for (yi, bi) in zip(y, bs))

        M = alpha.reduced_norm() // N
        if M <= B and M % 2 == 1:
            return M, alpha

    raise RuntimeError("random_bounded_elt failed")

def sample_quat_response(I_sk_chl, I_com, prng_res, prng_def, quat_ctx, ec_ctx):
    r"""
    [Algorithm 3.18] SampleQuaternionResponse.

    Sample the quaternion part of the response.

    INPUT:

    - ``I_sk_chl`` -- an ideal
    - ``I_com`` -- an ideal of norm coprime to ``I_sk_chl``
    - ``prng_res`` -- the response pseudorandom generator
    - ``prng_def`` -- the default-domain pseudorandom generator
    - ``quat_ctx`` -- the quaternion context
    - ``ec_ctx`` -- the elliptic curve context

    OUTPUT:

    A tuple ``(I_sk_chl1, alpha_rsp, alpha_sk_chl, q_rsp)`` with
    ``I_sk_chl1`` an ideal equivalent to ``I_sk_chl`` of prime norm and
    `q_{rsp} = \mathrm{nrd}(\alpha_{rsp}) / (\mathrm{nrd}(I_{com})
    \mathrm{nrd}(I'_{sk,chl}))`.
    """
    I_sk_chl = qt_inert.ideal_to_Nxy(I_sk_chl)

    I_sk_chl1, alpha_sk_chl = equivalent_coprime_ideal(
            I_sk_chl, 0, prng_def, quat_ctx)
    alpha_sk_chl = alpha_sk_chl.conjugate()

    if gcd(I_sk_chl1.norm(), I_com.norm()) != 1:
           raise RuntimeError("Commitment degree not coprime to sk_chall")

    q_rsp, alpha_rsp = random_bounded_elt(
        I_sk_chl1.conjugate()*I_com, 2**(ec_ctx.e_rsp-1), prng_res, quat_ctx)

    return I_sk_chl1, alpha_rsp, alpha_sk_chl, q_rsp

