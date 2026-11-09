from sage.all import (
    gcd, Integer, inverse_mod, log, Matrix, QQ, vector, xgcd, ZZ
)

def lat_to_basis(Lambda):
    r"""
    [Algorithm 4.139] LatToBasis.
    """
    return Lambda[:2], Lambda[2:]

def basis_to_lat(v, w):
    r"""
    [Algorithm 4.140] BasisToLat.
    """
    return v + w

def inner_product(v, w):
    r"""
    [Algorithm 4.136] InnerProduct.
    """
    return v[0]*w[0] + v[1]*w[1]

def vec_norm(v):
    r"""
    [Algorithm 4.135] Norm.
    """
    return inner_product(v, v)

def compute_2d_lattice(quat_ctx, N, alpha):
    r"""
    [Algorithm 4.141] Compute2DLattice.

    Return the kernel lattice of `(s, t) \mapsto \mathrm{tr}((s - it)\alpha)
    \bmod 4N`.

    INPUT:

    - ``quat_ctx`` -- the quaternion context
    - ``N`` -- a positive integer
    - ``alpha`` -- a quaternion element with
      `\gcd(\mathrm{tr}(\alpha), \mathrm{tr}(i\alpha), 2N) = 1`

    OUTPUT:

    A tuple ``(Lambda, v1)`` with ``Lambda`` the kernel lattice and ``v1``
    a 2D vector mapped to `1 \bmod 4N`.
    """
    a = alpha.reduced_trace()
    b = -(quat_ctx.i*alpha).reduced_trace()
    if gcd([a, b, 2*N]) != 1:
        raise ValueError("gcd(tr(alpha), tr(i*alpha), 2N) is not 1")
    N_four = 4*N
    g_a = gcd(a, N_four)
    a_m = a//g_a

    modulus = (N_four)//g_a
    c = inverse_mod(a_m, modulus)
    c = (-b * c) % modulus
    Lambda = [c, g_a, modulus, 0]
    k, u, v = xgcd(a, b)
    q, w, _ = xgcd(k, N_four)

    v1 = [(w * u) % N_four, (w * v) % N_four]

    return Lambda, v1

def one_step_reduction(v0, v1):
    r"""
    [Algorithm 4.144] OneStepReduction.

    INPUT:

    - ``v0``, ``v1`` -- 2D vectors with `\|v_0\| \le \|v_1\|`

    OUTPUT:

    A 2D vector `v_1 - q v_0` of norm at most `\|v_1\|`.
    """
    n0 = vec_norm(v0)
    if n0 > vec_norm(v1):
        raise ValueError("v0 is longer than v1")

    m = inner_product(v0, v1)

    q = round(m/n0)

    return list(vector(v1) - q*vector(v0))

def lagrange_reduction(Lambda):
    r"""
    Lagrange reduction of a 2D lattice. Together with
    :func:`compute_2d_lattice`, replaces [Algorithm 4.142]
    ComputeLatticeAndReduce.

    INPUT:

    - ``Lambda`` -- a 2D lattice

    OUTPUT:

    A reduced basis of ``Lambda``, as a 2D lattice.
    """

    v0, v1 = lat_to_basis(Lambda)
    if vec_norm(v0) > vec_norm(v1):
        v0, v1 = v1, v0

    found = False
    while not found:
        v1 = one_step_reduction(v0, v1)
        if vec_norm(v0) > vec_norm(v1):
            v0, v1 = v1, v0
        else:
            found = True

    return basis_to_lat(v0, v1)

def divide_minima(quat_ctx, N, beta_I, a, b, v):
    r"""
    [Algorithm 4.147] DivideMinima.

    INPUT:

    - ``quat_ctx`` -- the quaternion context
    - ``N`` -- a positive integer
    - ``beta_I`` -- an element of `O_0`
    - ``a``, ``b`` -- bounds with `0 < a < b`
    - ``v`` -- the smallest vector in the kernel lattice of ``beta_I``

    OUTPUT:

    A tuple ``(gamma_I, a_m, b_m, w)`` with ``gamma_I`` in
    `O_0 \langle \beta_I, N \rangle`, new bounds ``a_m < b_m`` and ``w`` the
    smallest vector in the kernel lattice of ``gamma_I``.
    """
    s0, t0 = v
    ell_I = gcd(s0, t0)

    B = log(N).n()**4
    g = log(B, 5).ceil()

    M = gcd(5**g, ell_I)
    h = Integer(M).valuation(5)

    s0_m = s0//M
    t0_m = t0//M

    if (2*s0_m - t0_m) % 5 != 0:
        tmp = (2 + quat_ctx.i)**h * (s0_m + quat_ctx.i*t0_m)
        s0_m, t0_m = list(tmp)[:2]
        gamma_I = (2 + quat_ctx.i)**h * beta_I
    else:
        tmp = (2 - quat_ctx.i)**h * (s0_m + quat_ctx.i*t0_m)
        s0_m, t0_m = list(tmp)[:2]
        gamma_I = (2 - quat_ctx.i)**h * beta_I

    a_m = M * a
    b_m = M * b
    return gamma_I, a_m, b_m, [s0_m, t0_m]

def reduce_given_minima(Lambda, w):
    r"""
    [Algorithm 4.145] ReduceGivenMinima.

    INPUT:

    - ``Lambda`` -- a 2D lattice
    - ``w`` -- a vector of smallest norm in ``Lambda``

    OUTPUT:

    A reduced basis of ``Lambda`` with first vector ``w``, as a 2D lattice.
    """

    if vector(ZZ, w) not in Matrix(ZZ, 2, 2, Lambda).row_module():
        raise ValueError("w is not in Lambda")
    a1, a2 = list(vector(w)*Matrix(QQ, 2, 2, Lambda)**(-1))
    gcd, a4, a3 = xgcd(a1, a2)

    if gcd != 1:
        raise ValueError("w is not primitive in Lambda")

    a3 = -a3

    v2 = list(vector([a3, a4])*Matrix(ZZ, 2, 2, Lambda))
    v2 = one_step_reduction(w, v2)

    return basis_to_lat(w, v2)

def compute_small_solutions(N, Lambda, lam, v, d):
    r"""
    [Algorithm 4.143] ComputeSmallSolutions.

    Find a small `(s, t)` with `\lambda(as + bt) \equiv d \bmod 4N`, using
    Babai rounding.

    INPUT:

    - ``N`` -- a positive integer
    - ``Lambda`` -- the kernel lattice of `(s, t) \mapsto as + bt \bmod 4N`
    - ``lam`` -- an integer coprime to `2N`
    - ``v`` -- a 2D vector with `as_1 + bt_1 \equiv 1 \bmod 4N`
    - ``d`` -- an integer

    OUTPUT:

    A tuple ``(s, t)`` of integers.
    """
    N_four = 4*N
    mu = inverse_mod(lam, N_four)

    s1, t1 = v
    vt = [(-mu*d*s1) % (N_four), (-mu*d*t1) % (N_four)]
    tmp = list(vector(vt)*Matrix(QQ, 2, 2, Lambda)**(-1))
    tmp = [round(c) for c in tmp]
    vc = vector(tmp) * Matrix(ZZ, 2, 2, Lambda)

    s, t = list(vc - vector(vt))
    return s, t
