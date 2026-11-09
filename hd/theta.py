def dim_two_dbl(P, U, V):
    r"""
    [Algorithm 4.114] Dim2DBL.

    INPUT:

    - ``P`` -- a :class:`DimTwoPoint`
    - ``U``, ``V`` -- the constants from :func:`compute_constants`

    OUTPUT:

    The point `[2]P`.
    """
    R = P.squaring()
    R = R.hadamard()
    R = R.squaring()
    R = R.pointwise_multiplication(V)
    R = R.hadamard()
    R = R.pointwise_multiplication(U)

    return R


def dim_two_dbl_iter(P, U, V, n):
    r"""
    Iterate :func:`dim_two_dbl` ``n`` times.

    INPUT:

    - ``P`` -- a :class:`DimTwoPoint`
    - ``U``, ``V`` -- the constants from :func:`compute_constants`
    - ``n`` -- a non-negative integer

    OUTPUT:

    The point `[2^n]P`.
    """
    R = P
    for _ in range(n):
        R = dim_two_dbl(R, U, V)
    return R


def compute_constants(A):
    r"""
    [Algorithm 4.112] ComputeConstants.

    INPUT:

    - ``A`` -- a :class:`DimTwoStructure`

    OUTPUT:

    A tuple ``(U, V)`` with `U = \theta_A(0_A)^{-1}` and
    `V = H(\theta_A(0_A)^2)^{-1}`.
    """
    OA = A.null_point()
    U = OA.inverse()
    V = OA.squaring()
    V = V.hadamard()
    V = V.inverse()
    return U, V


def recover_dual_null(U, x1, x2, y1, y2):
    r"""
    Recover the theta structure from ``U`` and the data of
    :func:`generic_codomain`. First part of [Algorithm 4.113].

    INPUT:

    - ``U`` -- the point `\theta_A(0_A)^{-1}`
    - ``x1``, ``x2``, ``y1``, ``y2`` -- the constants from
      :func:`generic_codomain`

    OUTPUT:

    The :class:`DimTwoStructure` with null point `\theta_A(0_A)`.
    """
    _, _, c_inv, d_inv = U.dim_two_coordinates()

    u = x2 * y2
    a_prime = u * x1
    b_prime = u * y1
    c_prime = d_inv
    d_prime = c_inv

    OA = DimTwoPoint([a_prime, b_prime, c_prime, d_prime])

    return DimTwoStructure(OA)


def compute_constants_from_iso_data(U, x1, x2, y1, y2):
    r"""
    [Algorithm 4.113] ComputeConstantsFromIsoData.

    INPUT:

    - ``U`` -- the point `\theta_A(0_A)^{-1}`
    - ``x1``, ``x2``, ``y1``, ``y2`` -- the constants from
      :func:`generic_codomain`

    OUTPUT:

    A tuple ``(A, V)`` with ``A`` the :class:`DimTwoStructure` and
    `V = H(\theta_A(0_A)^2)^{-1}`.
    """
    A = recover_dual_null(U, x1, x2, y1, y2)
    OA = A.null_point()

    V = OA.squaring()
    V = V.hadamard()
    V = V.inverse()

    return A, V


def verify_kernel(P, U):
    r"""
    [Algorithm 4.115] VerifyKernel.

    INPUT:

    - ``P`` -- the point `H(\theta_A(K_1)^2)`
    - ``U`` -- the point `\tilde\theta_B(0_B)^{-1}`

    OUTPUT:

    ``True`` if the kernel is maximal isotropic, ``False`` otherwise.
    """
    _, _, zP, tP = P.dim_two_coordinates()
    _, _, zU, tU = U.dim_two_coordinates()

    t1 = zP * zU
    t2 = tP * tU
    if t1 != t2:
        return False
    return True


def generic_codomain(K1, K2, dual_codomain, check_kernel):
    r"""
    [Algorithm 4.116] GenericCodomain.

    Compute the codomain of the `(2, 2)`-isogeny with kernel
    `\langle [4]K_1, [4]K_2 \rangle`.

    INPUT:

    - ``K1``, ``K2`` -- points of type :class:`DimTwoPoint` of order 8
    - ``dual_codomain`` -- if ``True``, apply the Hadamard transform first
    - ``check_kernel`` -- if ``True``, verify the kernel

    OUTPUT:

    A tuple ``(U, x1, x2, y1, y2)`` with `U = \tilde\theta_B(0_B)^{-1}` and
    the constants for :func:`compute_constants_from_iso_data`.
    """
    if dual_codomain:
        P1 = K1.hadamard()
        P2 = K2.hadamard()
    else:
        P1 = K1
        P2 = K2

    P1 = P1.squaring()
    P2 = P2.squaring()
    P1 = P1.hadamard()
    P2 = P2.hadamard()

    x1, y1, _, _ = P1.dim_two_coordinates()
    x2, y2, z2, t2 = P2.dim_two_coordinates()

    if check_kernel:
        if 0 in [x2, y2, z2, t2, x1, y1]:
            raise ValueError("generic codomain failed: unexpected splitting occurred")

    u1 = y1 * t2
    u2 = x1 * z2
    A_inv = u1 * z2
    B_inv = u2 * t2
    C_inv = u1 * x2
    D_inv = u2 * y2

    U = DimTwoPoint([A_inv, B_inv, C_inv, D_inv])

    if check_kernel:
        if not verify_kernel(P1, U):
            raise ValueError("generic codomain failed: kernel is invalid")

    return U, x1, x2, y1, y2


def generic_eval(R, U, dual_codomain):
    r"""
    [Algorithm 4.117] GenericEval, on a single point.

    INPUT:

    - ``R`` -- a :class:`DimTwoPoint`
    - ``U`` -- the output of :func:`generic_codomain`
    - ``dual_codomain`` -- as in :func:`generic_codomain`

    OUTPUT:

    The image of ``R``, as a :class:`DimTwoPoint`.
    """
    P = R
    if dual_codomain:
        P = P.hadamard()

    P = P.squaring()
    P = P.hadamard()
    P = P.pointwise_multiplication(U)

    return P


class DimTwoPoint:
    r"""
    Implements the type Dim2Point [Definition 4.5.5].

    A theta point `(x : y : z : t) \in \mathbb{P}^3(\mathbb{F}_{p^2})`.

    INPUT:

    - ``coords`` -- the list ``[x, y, z, t]``
    """

    def __init__(self, coords):
        self._coords = coords

    def dim_two_coordinates(self):
        r"""
        Return the coordinates ``[x, y, z, t]`` [Definition 4.5.5].
        """
        return self._coords

    def pointwise_multiplication(self, other):
        r"""
        [Algorithm 4.108] PointwiseMult.

        Return the coordinate-wise product of ``self`` and ``other``.
        """
        x1, y1, z1, t1 = self._coords
        x2, y2, z2, t2 = other._coords

        X = x1 * x2
        Y = y1 * y2
        Z = z1 * z2
        T = t1 * t2

        return DimTwoPoint([X, Y, Z, T])

    def squaring(self):
        r"""
        [Algorithm 4.109] Squaring.

        Return the coordinate-wise square of ``self``.
        """
        x, y, z, t = self._coords

        X = x * x
        Y = y * y
        Z = z * z
        T = t * t

        return DimTwoPoint([X, Y, Z, T])

    def inverse(self):
        r"""
        [Algorithm 4.110] Inverse.

        Return the coordinate-wise projective inverse of ``self``.
        """
        x, y, z, t = self._coords

        t1 = x * y
        t2 = z * t

        X = t2 * y
        Y = t2 * x
        Z = t1 * t
        T = t1 * z

        return DimTwoPoint([X, Y, Z, T])

    def hadamard(self):
        r"""
        [Algorithm 4.111] Hadamard.

        Return the Hadamard transform of ``self``.
        """
        x, y, z, t = self._coords

        t1 = x + y
        t2 = x - y
        t3 = z + t
        t4 = z - t

        X = t1 + t3
        Y = t2 + t4
        Z = t1 - t3
        T = t2 - t4

        return DimTwoPoint([X, Y, Z, T])

    def __repr__(self):
        return f"{list(self._coords)}"

    def __str__(self):
        return f"ThetaPoint: {list(self._coords)}"

    def __eq__(self, other):
        r"""
        Projective equality of two points.
        """
        if not isinstance(other, DimTwoPoint):
            return False

        a1, b1, c1, d1 = self._coords
        a2, b2, c2, d2 = other._coords

        if d1 != 0 or d2 != 0:
            return all([a1 * d2 == a2 * d1, b1 * d2 == b2 * d1, c1 * d2 == c2 * d1])
        elif c1 != 0 or c2 != 0:
            return all([a1 * c2 == a2 * c1, b1 * c2 == b2 * c1])
        elif b1 != 0 or b2 != 0:
            return a1 * b2 == a2 * b1
        else:
            return True


class DimTwoStructure:
    r"""
    Implements the type Dim2Structure [Definition 4.5.4].

    A theta null point `(a : b : c : d) \in \mathbb{P}^3(\mathbb{F}_{p^2})`.

    INPUT:

    - ``null_point`` -- the null point, as a :class:`DimTwoPoint`
    """

    def __init__(self, null_point):
        self._null_point = null_point

    def null_point(self):
        r"""
        Return the null point, as a :class:`DimTwoPoint`.
        """
        return self._null_point
