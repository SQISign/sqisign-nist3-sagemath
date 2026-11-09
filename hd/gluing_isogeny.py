from sage.all import Matrix

from hd.couple_point import CoupleECPoint
from hd.theta import DimTwoPoint, DimTwoStructure


def verify_two_torsion(P, Q):
    r"""
    Check that ``P`` and ``Q`` generate a `(2, 2)`-subgroup of
    `E_1 \times E_2`.

    INPUT:

    - ``P``, ``Q`` -- points of type :class:`CoupleECPoint`

    OUTPUT:

    ``True`` if all components have order 2 and ``P``, ``Q`` are independent
    componentwise, ``False`` otherwise.
    """
    P1, P2 = P.points()
    Q1, Q2 = Q.points()

    if not all(R.has_order_2e(1) for R in (P1, P2, Q1, Q2)):
        return False

    # Now ensure independence
    if P1 == Q1 or P2 == Q2:
        return False

    return True


def co_sym_elt(P, Q):
    r"""
    [Algorithm 4.118] CoSymElem.

    INPUT:

    - ``P``, ``Q`` -- a basis of the 4-torsion of a Montgomery curve

    OUTPUT:

    A tuple ``(gP, gQ, gPQ, lam)`` with the `2 \times 2` matrices associated
    to `P`, `Q`, `P + Q` and the projective factor ``lam``.
    """
    xP, zP = P.xz()
    xQ, zQ = Q.xz()

    i = xP.parent().gen()

    sP = xP + zP
    sQ = xQ + zQ
    tP = xP - zP
    tQ = xQ - zQ

    tP = tP * sQ
    tQ = tQ * sP
    tQ = i * tQ

    if xP == zP:
        x2 = xQ**2
        z2 = zQ**2

        lam = x2 - z2
        r1 = x2 + z2
        r2 = xQ * zQ
        r2 = r2 + r2

        s1 = i * r1
        s2 = i * r2

        gT1 = Matrix([[0, lam], [lam, 0]])
        gT2 = Matrix([[r1, -r2], [r2, -r1]])
        gT12 = Matrix([[s2, -s1], [s1, -s2]])

    elif xP == -zP:
        x2 = xQ**2
        z2 = zQ**2

        lam = x2 - z2
        r1 = x2 + z2
        r2 = xQ * zQ
        r2 = r2 + r2

        s1 = i * r1
        s2 = i * r2

        gT1 = Matrix([[0, -lam], [-lam, 0]])
        gT2 = Matrix([[r1, -r2], [r2, -r1]])
        gT12 = Matrix([[-s2, s1], [-s1, s2]])

    elif xQ == zQ:
        x2 = xP**2
        z2 = zP**2

        lam = x2 - z2
        r1 = x2 + z2
        r2 = xP * zP
        r2 = r2 + r2

        s1 = i * r1
        s2 = i * r2

        gT1 = Matrix([[r1, -r2], [r2, -r1]])
        gT2 = Matrix([[0, lam], [lam, 0]])
        gT12 = Matrix([[-s2, s1], [-s1, s2]])

    elif xQ == -zQ:
        x2 = xP**2
        z2 = zP**2

        lam = x2 - z2
        r1 = x2 + z2
        r2 = xP * zP
        r2 = r2 + r2

        s1 = i * r1
        s2 = i * r2

        gT1 = Matrix([[r1, -r2], [r2, -r1]])
        gT2 = Matrix([[0, -lam], [-lam, 0]])
        gT12 = Matrix([[s2, -s1], [s1, -s2]])

    elif tQ == -tP:
        x2 = xP**2
        z2 = zP**2

        lam = x2 - z2
        r1 = x2 + z2
        r2 = xP * zP
        r2 = r2 + r2

        s1 = i * r1
        s2 = i * r2

        gT1 = Matrix([[r1, -r2], [r2, -r1]])
        gT2 = Matrix([[-s2, s1], [-s1, s2]])
        gT12 = Matrix([[0, -lam], [-lam, 0]])

    elif tQ == tP:
        x2 = xP**2
        z2 = zP**2

        lam = x2 - z2
        r1 = x2 + z2
        r2 = xP * zP
        r2 = r2 + r2

        s1 = i * r1
        s2 = i * r2

        gT1 = Matrix([[r1, -r2], [r2, -r1]])
        gT2 = Matrix([[s2, -s1], [s1, -s2]])
        gT12 = Matrix([[0, lam], [lam, 0]])
    else:
        raise ValueError(
            "CoSymElt failed: input points do not form a basis of the 4-torsion of a Montgomery curve."
        )

    return gT1, gT2, gT12, lam


def theta_change_of_coordinates(K1, K2):
    r"""
    [Algorithm 4.119] ThetaChangeOfCoordinates.

    INPUT:

    - ``K1``, ``K2`` -- points of type :class:`CoupleECPoint` of order 4,
      such that `\langle [2]K_1, [2]K_2 \rangle` is the gluing kernel

    OUTPUT:

    The `4 \times 4` change of coordinates matrix.
    """
    K11, K12 = K1.points()
    K21, K22 = K2.points()

    b1_g1, b1_g2, b1_g3, lam_1 = co_sym_elt(K11, K21)
    b2_g1, b2_g2, b2_g3, lam_2 = co_sym_elt(K12, K22)
    lam_12 = lam_1 * lam_2

    M00 = (
        lam_12
        + b1_g1[0][0] * b2_g1[0][0]
        + b1_g2[0][0] * b2_g2[0][0]
        - b1_g3[0][0] * b2_g3[0][0]
    )
    M01 = (
        b1_g1[0][0] * b2_g1[0][1]
        + b1_g2[0][0] * b2_g2[0][1]
        - b1_g3[0][0] * b2_g3[0][1]
    )
    M02 = (
        b1_g1[0][1] * b2_g1[0][0]
        + b1_g2[0][1] * b2_g2[0][0]
        - b1_g3[0][1] * b2_g3[0][0]
    )
    M03 = (
        b1_g1[0][1] * b2_g1[0][1]
        + b1_g2[0][1] * b2_g2[0][1]
        - b1_g3[0][1] * b2_g3[0][1]
    )

    M10 = M00 * b2_g2[0][0] + M01 * b2_g2[1][0]
    M11 = M00 * b2_g2[0][1] - M01 * b2_g2[0][0]
    M12 = M02 * b2_g2[0][0] + M03 * b2_g2[1][0]
    M13 = M02 * b2_g2[0][1] - M03 * b2_g2[0][0]

    M20 = M00 * b1_g1[0][0] + M02 * b1_g1[1][0]
    M21 = M01 * b1_g1[0][0] + M03 * b1_g1[1][0]
    M22 = M00 * b1_g1[0][1] - M02 * b1_g1[0][0]
    M23 = M01 * b1_g1[0][1] - M03 * b1_g1[0][0]

    M30 = M10 * b1_g1[0][0] + M12 * b1_g1[1][0]
    M31 = M11 * b1_g1[0][0] + M13 * b1_g1[1][0]
    M32 = M10 * b1_g1[0][1] - M12 * b1_g1[0][0]
    M33 = M11 * b1_g1[0][1] - M13 * b1_g1[0][0]

    M00 = lam_12 * M00
    M01 = lam_12 * M01
    M02 = lam_12 * M02
    M03 = lam_12 * M03

    M10 = lam_1 * M10
    M11 = lam_1 * M11
    M12 = lam_1 * M12
    M13 = lam_1 * M13

    M20 = lam_2 * M20
    M21 = lam_2 * M21
    M22 = lam_2 * M22
    M23 = lam_2 * M23

    M = Matrix(4, 4, [
        M00, M01, M02, M03,
        M10, M11, M12, M13,
        M20, M21, M22, M23,
        M30, M31, M32, M33,
    ])

    return M


def apply_matrix(M, P):
    r"""
    [Algorithm 4.120] ApplyMatrix.

    INPUT:

    - ``M`` -- a `4 \times 4` matrix
    - ``P`` -- a :class:`DimTwoPoint`

    OUTPUT:

    The :class:`DimTwoPoint` `M P^T`.
    """
    xP, yP, zP, tP = P.dim_two_coordinates()
    xQ = M[0][0] * xP + M[0][1] * yP + M[0][2] * zP
    yQ = M[1][0] * xP + M[1][1] * yP + M[1][2] * zP
    zQ = M[2][0] * xP + M[2][1] * yP + M[2][2] * zP
    tQ = M[3][0] * xP + M[3][1] * yP + M[3][2] * zP
    if tP != 0:
        xQ = xQ + M[0][3] * tP
        yQ = yQ + M[1][3] * tP
        zQ = zQ + M[2][3] * tP
        tQ = tQ + M[3][3] * tP
    return DimTwoPoint([xQ, yQ, zQ, tQ])


def product_to_theta(M, P):
    r"""
    [Algorithm 4.121] ProductToTheta.

    INPUT:

    - ``M`` -- the change of coordinates matrix
    - ``P`` -- a point of type :class:`CoupleECPoint`

    OUTPUT:

    The theta coordinates of ``P``, as a :class:`DimTwoPoint`.
    """
    P1, P2 = P.points()
    x1, z1 = P1.xz()
    x2, z2 = P2.xz()

    xQ = x1 * x2
    yQ = x1 * z2
    zQ = z1 * x2
    tQ = z1 * z2

    Q = DimTwoPoint([xQ, yQ, zQ, tQ])
    Q = apply_matrix(M, Q)

    return Q


def gluing_codomain(K1, K2, dual_codomain, verify):
    r"""
    [Algorithm 4.122] GluingCodomain.

    Compute the codomain of the gluing isogeny with kernel
    `\langle [4]K_1, [4]K_2 \rangle`.

    INPUT:

    - ``K1``, ``K2`` -- points of type :class:`CoupleECPoint` of order 8
    - ``dual_codomain`` -- if ``False``, apply the Hadamard transform to the
      codomain null point
    - ``verify`` -- if ``True``, check that `[4]K_1, [4]K_2` generate a
      `(2, 2)`-subgroup

    OUTPUT:

    A tuple ``(B, U, I1, I2, M)`` with ``B`` the codomain, ``U`` the inverse
    dual theta null point, ``I1``, ``I2`` the inverse images of ``K1``,
    ``K2`` and ``M`` the change of coordinates matrix.
    """
    K1_prime = K1.dbl()
    K2_prime = K2.dbl()

    if verify:
        K1_2 = K1_prime.dbl()
        K2_2 = K2_prime.dbl()
        if not verify_two_torsion(K1_2, K2_2):
            raise ValueError("not two torsion")

    M = theta_change_of_coordinates(K1_prime, K2_prime)

    P1 = product_to_theta(M, K1)
    P1 = P1.squaring()
    P1 = P1.hadamard()
    x1, y1, z1, t1 = P1.dim_two_coordinates()

    P2 = product_to_theta(M, K2)
    P2 = P2.squaring()
    P2 = P2.hadamard()
    x2, y2, z2, t2 = P2.dim_two_coordinates()

    if t1 != 0 or t2 != 0:
        raise ValueError("GluingCodomain failed: kernel is not isotropic")

    A = x1 * x2
    B = y1 * x2
    C = x1 * z2

    A_inv = y1 * z2
    B_inv = C
    C_inv = B

    x_inv = z1 * C_inv
    y_inv = x1 * A_inv
    z_inv = y2 * B_inv
    t_inv = x2 * A_inv

    zero = A.parent().zero()

    OB = DimTwoPoint([A, B, C, zero])
    if not dual_codomain:
        OB = OB.hadamard()
    B = DimTwoStructure(OB)

    U = DimTwoPoint([A_inv, B_inv, C_inv, zero])
    I1 = DimTwoPoint([x_inv, x_inv, y_inv, y_inv])
    I2 = DimTwoPoint([z_inv, t_inv, z_inv, t_inv])

    return B, U, I1, I2, M


def gluing_eval(P1, P2, I, M, dual_codomain):
    r"""
    [Algorithm 4.123] GluingEval.

    Evaluate the gluing isogeny at `(P_1, P_2)`, given in barycentric
    coordinates with respect to `W = (W_1, W_2)`.

    INPUT:

    - ``P1``, ``P2`` -- barycentric coordinates ``(u, v, w)`` of
      `P_1 \pm W_1` and `P_2 \pm W_2`
    - ``I`` -- the inverse dual theta coordinates of the image of `W`
    - ``M`` -- the change of coordinates matrix
    - ``dual_codomain`` -- if ``True``, skip the final Hadamard transform

    OUTPUT:

    The image of `(P_1, P_2)`, in dual theta coordinates if
    ``dual_codomain`` is ``True``.
    """
    u1, v1, w1 = P1
    u2, v2, w2 = P2

    xU = u1 * u2
    tV = v1 * v2
    xU = xU + tV
    yU = u1 * w2
    zU = w1 * u2
    tU = w1 * w2

    xV = u1 + v1
    yV = u2 + v2
    xV = xV * yV
    xV = xV - xU
    yV = v1 * w2
    zV = w1 * v2

    zero = u1.parent().zero()
    U = DimTwoPoint([xU, yU, zU, tU])
    V = DimTwoPoint([xV, yV, zV, zero])

    U = apply_matrix(M, U)
    V = apply_matrix(M, V)

    U = U.squaring()
    V = V.squaring()

    xU, yU, zU, tU = U.dim_two_coordinates()
    xV, yV, zV, tV = V.dim_two_coordinates()
    x = xU - xV
    y = yU - yV
    z = zU - zV
    t = tU - tV

    Q = DimTwoPoint([x, y, z, t])
    Q = Q.hadamard()
    Q = Q.pointwise_multiplication(I)

    if not dual_codomain:
        Q = Q.hadamard()
    return Q


def gluing_eval_special(P, U, M):
    r"""
    [Algorithm 4.124] GluingEvalSpecial.

    INPUT:

    - ``P`` -- a point of type :class:`CoupleECPoint` of the form
      `(P_1, 0)` or `(0, P_2)`
    - ``U`` -- the inverse dual theta null point of the codomain
    - ``M`` -- the change of coordinates matrix

    OUTPUT:

    The image of ``P``, as a :class:`DimTwoPoint`.
    """
    Q = product_to_theta(M, P)
    Q = Q.squaring()
    Q = Q.hadamard()
    xQ, yQ, zQ, tQ = Q.dim_two_coordinates()
    if tQ != 0:
        raise ValueError("GluingEvalSpecial: last coordinate should be zero")

    alpha_inv, beta_inv, gamma_inv, _ = U.dim_two_coordinates()

    xQ = xQ * alpha_inv
    yQ = yQ * beta_inv
    zQ = zQ * gamma_inv

    zero = xQ.parent().zero()
    Q = DimTwoPoint([xQ, yQ, zQ, zero])
    return Q.hadamard()


def first_isogeny(E1, E2, T1, T2, e, pts, verify):
    r"""
    [Algorithm 4.131] FirstIsogeny.

    Compute the first (gluing) step `\Phi_1 : E_1 \times E_2 \to B` of a
    `(2, 2)`-isogeny chain of length ``e``.

    INPUT:

    - ``E1``, ``E2`` -- the domain curves
    - ``T1``, ``T2`` -- points of type :class:`CoupleECPoint` generating the
      kernel of the chain
    - ``e`` -- the length of the chain
    - ``pts`` -- a list of :class:`CoupleECPoint` of the form `(P_1, 0)` or
      `(0, P_2)`
    - ``verify`` -- if ``True``, check the gluing kernel

    OUTPUT:

    A tuple ``(B, orders, theta_q1, theta_q2, image_pts, k)`` with ``B`` the
    codomain, ``theta_q1``, ``theta_q2`` the images of the strategy points
    `[2^{e - r}]T_1`, `[2^{e - r}]T_2` for `r` in ``orders``, ``image_pts``
    the images of ``pts`` and ``k`` the number of strategy points.
    """
    T11, T12 = T1.points()
    T21, T22 = T2.points()

    space = e.bit_length() + 1
    orders = [0] * space

    ker_11_even = [None] * space
    ker_12_even = [None] * space
    ker_21_even = [None] * space
    ker_22_even = [None] * space
    ker_11_odd = [None] * space
    ker_12_odd = [None] * space
    ker_21_odd = [None] * space
    ker_22_odd = [None] * space

    theta_q1 = [None] * space
    theta_q2 = [None] * space

    orders[0] = e
    ker_11_even[0] = T11
    ker_12_even[0] = T12
    ker_21_even[0] = T21
    ker_22_even[0] = T22
    ker_11_odd[0] = E1.zero()
    ker_12_odd[0] = E2.zero()
    ker_21_odd[0] = E1.zero()
    ker_22_odd[0] = E2.zero()

    k = 0
    while orders[k] != 1:
        if orders[k] >= 16:
            n = orders[k] // 2
        else:
            n = orders[k] - 1

        ker_11_even[k + 1], ker_11_odd[k + 1] = ker_11_even[k].ec_gluing_ladder(
            ker_11_odd[k], T11, n
        )
        ker_12_even[k + 1], ker_12_odd[k + 1] = ker_12_even[k].ec_gluing_ladder(
            ker_12_odd[k], T12, n
        )
        ker_21_even[k + 1], ker_21_odd[k + 1] = ker_21_even[k].ec_gluing_ladder(
            ker_21_odd[k], T21, n
        )
        ker_22_even[k + 1], ker_22_odd[k + 1] = ker_22_even[k].ec_gluing_ladder(
            ker_22_odd[k], T22, n
        )

        orders[k + 1] = orders[k] - n

        k = k + 1

    if not all(R.has_order_2e(3) for R in (ker_11_even[k], ker_12_even[k],
                                            ker_21_even[k], ker_22_even[k])):
        raise ValueError("FirstIsogeny failed: kernel points have wrong order")

    K1 = CoupleECPoint(ker_11_even[k], ker_12_even[k])
    K2 = CoupleECPoint(ker_21_even[k], ker_22_even[k])

    B, U, I1, I2, M = gluing_codomain(K1, K2, False, verify)

    K11, K12 = K1.points()
    K21, K22 = K2.points()

    T11_K11 = T11.ec_bary_components(K11, ker_11_odd[k])
    T12_K12 = T12.ec_bary_components(K12, ker_12_odd[k])
    T21_K21 = T21.ec_bary_components(K21, ker_21_odd[k])
    T22_K22 = T22.ec_bary_components(K22, ker_22_odd[k])

    I1 = gluing_eval(T11_K11, T12_K12, I1, M, True)
    I2 = gluing_eval(T21_K21, T22_K22, I2, M, True)

    theta_q1[0] = I1.hadamard()
    theta_q2[0] = I2.hadamard()
    orders[0] -= 1

    I1_inv = I1.inverse()
    I2_inv = I2.inverse()

    for j in range(1, k):
        K11_T11 = ker_11_even[j].ec_bary_components(T11, ker_11_odd[j])
        K12_T12 = ker_12_even[j].ec_bary_components(T12, ker_12_odd[j])
        K21_T21 = ker_21_even[j].ec_bary_components(T21, ker_21_odd[j])
        K22_T22 = ker_22_even[j].ec_bary_components(T22, ker_22_odd[j])

        theta_q1[j] = gluing_eval(K11_T11, K12_T12, I1_inv, M, False)
        theta_q2[j] = gluing_eval(K21_T21, K22_T22, I2_inv, M, False)
        orders[j] -= 1

    image_pts = [gluing_eval_special(P, U, M) for P in pts]
    k = k - 1

    return B, orders, theta_q1, theta_q2, image_pts, k
