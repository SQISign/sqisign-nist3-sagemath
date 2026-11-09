from sage.all import ceil, isqrt, is_prime, gcd, ZZ, Integers, kronecker


def coordinates_O0(alpha):
    r"""
    [Algorithm 4.43] CoordinatesO0.

    Return the coordinates of ``alpha`` in the `\ZZ`-basis
    `(1, i, (i+j)/2, (1+k)/2)` of the maximal order `O_0`.

    INPUT:

    - ``alpha`` -- a quaternion element of `O_0`

    OUTPUT:

    A list ``[a, b, c, d]`` of integers such that ``alpha = a + b*i + c*(i+j)/2
    + d*(1+k)/2``.
    """
    x1, x2, x3, x4 = list(alpha)

    a = x1 - x4
    b = x2 - x3
    c = 2 * x3
    d = 2 * x4

    coeffs = [a, b, c, d]

    if not all(r in ZZ for r in coeffs):
        raise ValueError("alpha is not in O0")
    return coeffs

def mod_O0(alpha, N):
    r"""
    [Algorithm 4.44] ModO0.

    Reduce the `O_0`-coordinates of ``alpha`` modulo ``N``.

    INPUT:

    - ``alpha`` -- a quaternion element of `O_0`
    - ``N`` -- a positive integer

    OUTPUT:

    A tuple ``(a, b, c, d)`` with entries in `[0, N)`, the coordinates of
    ``alpha`` in the basis `(1, i, (i+j)/2, (1+k)/2)` reduced modulo ``N``.
    """
    a, b, c, d = coordinates_O0(alpha)
    return a % N, b % N, c % N, d % N

def from_O0_coords(alpha, quat_ctx):
    r"""
    Return the quaternion element with the given `O_0`-coordinates.

    INPUT:

    - ``alpha`` -- a sequence ``(a, b, c, d)`` of integers
    - ``quat_ctx`` -- the quaternion context

    OUTPUT:

    The element ``a + b*i + c*(i+j)/2 + d*(1+k)/2``.
    """
    i, j, k = quat_ctx.B.gens()
    a, b, c, d = alpha
    return a + b*i + c*(i+j)/2 + d*(1+k)/2

def coordinates_B(alpha):
    r"""
    [Algorithm 4.41] CoordinatesB.

    Return the coordinates of ``alpha`` in the basis `(1, i, j, k)`, as
    integer numerators plus a common denominator.

    INPUT:

    - ``alpha`` -- a quaternion element

    OUTPUT:

    A tuple ``(a, b, c, d, e)`` of integers such that
    ``alpha = (a + b*i + c*j + d*k)/e``.
    """
    a, b, c, d = list(alpha)
    e = alpha.denominator()
    return a * e, b * e, c * e, d * e, e

def cornacchia(m):
    r"""
    [Algorithm 4.49] Cornacchia.

    Solve `x^2 + y^2 = m` for a prime ``m``.

    INPUT:

    - ``m`` -- a prime integer

    OUTPUT:

    A tuple ``(x, y)`` of integers with `x^2 + y^2 = m` and `0 \le x \le y`,
    or ``False`` if no solution exists.
    """
    if kronecker(-1, m) != 1:
        return False

    if m == 2:
        return 1, 1

    r = ZZ(Integers(m)(-1).sqrt())
    s = m

    r, s = s, r % s
    while True:
        if s**2 < m <= r**2:
            break
        r, s = s, r % s

    x, y = s, isqrt(m - s**2)
    if y < x:
        x, y = y, x

    if x**2 + y**2 == m:
        return x, y


def represent_integer(M, prng, quat_ctx):
    r"""
    [Algorithm 4.46] RepresentInteger.

    Find a primitive element of `O_0` of reduced norm ``M``.

    INPUT:

    - ``M`` -- an odd integer with `M > p`
    - ``prng`` -- the pseudorandom generator
    - ``quat_ctx`` -- the quaternion context

    OUTPUT:

    A quaternion element ``gamma`` of `O_0` with `\mathrm{nrd}(\gamma) = M`.
    """

    bound = ceil(M / quat_ctx.p)
    counter = 0

    while counter < bound:
        counter += 1
        m = isqrt(M / quat_ctx.p - 1)
        z = prng.sample_from_interval(1, m)

        m1 = isqrt((M - quat_ctx.p * z**2) / quat_ctx.p)
        t = prng.sample_from_interval(1, m1)

        M1 = ZZ(M - quat_ctx.p * (z**2 + t**2))

        if is_prime(M1):
            res = cornacchia(M1)
            if not res:
                continue
            x, y = res

            g = gcd([x, y, z, t])
            if g != 1:
                continue

            return quat_ctx.B([x, y, z, t])
    raise RuntimeError("represent_integer failed")

