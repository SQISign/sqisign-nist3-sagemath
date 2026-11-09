from ec.ec_curve import ECCurve
from hd.theta import DimTwoPoint
from hd.couple_point import CoupleECPoint
import crypto.utils as utils

def compute_mont_coeff(a, b):
    r"""
    [Algorithm 4.125] ComputeMontCoeff.

    INPUT:

    - ``a``, ``b`` -- a one-dimensional theta null point

    OUTPUT:

    A tuple ``(A, C)`` defining projectively the Montgomery coefficient.
    """
    a4 = a * a
    a4 = a4 * a4
    b4 = b * b
    b4 = b4 * b4

    A = a4 + b4
    A = A + A
    C = b4 - a4

    if C == 0:
        raise ValueError("ComputeMontCoeff failed with C=0")

    return (A, C)


def theta_to_montgomery(a, b):
    r"""
    [Algorithm 4.126] ThetaToMontgomery.

    INPUT:

    - ``a``, ``b`` -- a one-dimensional theta null point

    OUTPUT:

    A tuple ``(E, N)`` with ``E`` the curve with maximal Montgomery
    coefficient and ``N`` the `2 \times 2` matrix, as a list, converting theta
    coordinates into Montgomery `(x : z)` coordinates on ``E``.
    """
    i = a.parent().gen()

    a_prime = i * a
    b_prime = i * b

    a1 = a + b
    b1 = a - b
    a2 = a + b_prime
    b2 = a - b_prime

    A0, C0 = compute_mont_coeff(a, b)
    A1, C1 = compute_mont_coeff(a1, b1)
    A2, C2 = compute_mont_coeff(a2, b2)
    C0, C1, C2 = utils.batched_inversion(C0, C1, C2)

    A0 *= C0
    A1 *= C1
    A2 *= C2
    A3 = -A0
    A4 = -A1
    A5 = -A2

    Ai = [A0, A1, A2, A3, A4, A5]
    A_max = Ai[0]
    i_max = 0
    for i in range(1, 6):
        if utils.fp2_compare(A_max, Ai[i]):
            i_max = i
            A_max = Ai[i]

    if i_max == 0:
        N = [b, a, -b, a]
    elif i_max == 1:
        N = [a, -b, b, -a]
    elif i_max == 2:
        N = [a, b, b_prime, -a_prime]
    elif i_max == 3:
        N = [b, a, b, -a]
    elif i_max == 4:
        N = [a, -b, -b, a]
    elif i_max == 5:
        N = [a, b, -b_prime, a_prime]

    return ECCurve(A_max), N


def non_prod_to_prod(P, hadamard):
    r"""
    [Algorithm 4.127] NonProdToProd.

    INPUT:

    - ``P`` -- a :class:`DimTwoPoint`
    - ``hadamard`` -- a boolean

    OUTPUT:

    The point ``P`` in product theta coordinates.
    """
    x, y, z, t = P.dim_two_coordinates()

    if hadamard:
        x_prime = x + z
        y_prime = y + t
        z_prime = y - t
        t_prime = x - z
    else:
        x_prime = x + y
        y_prime = x - y
        z_prime = z - t
        t_prime = z + t

    return DimTwoPoint([x_prime, y_prime, z_prime, t_prime])


def dim_two_struct_to_ec_curves(A, hadamard):
    r"""
    [Algorithm 4.128] Dim2StructToECCurves.

    INPUT:

    - ``A`` -- a :class:`DimTwoStructure` isomorphic to a product `E_1 \times
      E_2`
    - ``hadamard`` -- a boolean

    OUTPUT:

    A tuple ``(E1, E2, N1, N2)`` with the curves and the matrices converting
    theta coordinates into Montgomery coordinates on each of them.
    """
    OA = A.null_point()
    OA_prime = non_prod_to_prod(OA, hadamard)
    a_prime, b_prime, c_prime, d_prime = OA_prime.dim_two_coordinates()
    if 0 in [a_prime, b_prime, c_prime, d_prime]:
        raise ValueError("zero found")

    t1 = a_prime * d_prime
    t2 = b_prime * c_prime

    if a_prime == 0 or b_prime == 0 or c_prime == 0 or t1 != t2:
        raise ValueError(
            "DimTwoStructToECCurves failed: input not isomorphic to a product of Montgomery curves"
        )

    E1, N1 = theta_to_montgomery(a_prime, b_prime)
    E2, N2 = theta_to_montgomery(a_prime, c_prime)

    return E1, E2, N1, N2


def dim_two_point_to_ec_points(P, E1, E2, N1, N2, hadamard):
    r"""
    [Algorithm 4.129] DimTwoPointToECPoints.

    INPUT:

    - ``P`` -- a :class:`DimTwoPoint`
    - ``E1``, ``E2`` -- the curves of the product
    - ``N1``, ``N2`` -- the conversion matrices for ``E1`` and ``E2``
    - ``hadamard`` -- a boolean

    OUTPUT:

    The point ``P`` as a :class:`CoupleECPoint` on `E_1 \times E_2`.
    """
    Q = non_prod_to_prod(P, hadamard)
    x, y, z, t = Q.dim_two_coordinates()

    if x != 0:
        th_10 = x
        th_11 = y
        th_20 = x
        th_21 = z
    elif y != 0:
        th_10 = x
        th_11 = y
        th_20 = y
        th_21 = t
    elif z != 0:
        th_10 = z
        th_11 = t
        th_20 = x
        th_21 = z
    else:
        raise ValueError("Dim2StructToECCurves fails: zero input")

    t1 = x * t
    t2 = y * z

    if t1 != t2:
        raise ValueError("Dim2StructToECCurves fails: input is not a product")

    N1_00, N1_01, N1_10, N1_11 = N1
    N2_00, N2_01, N2_10, N2_11 = N2

    x1 = N1_00 * th_10 + N1_01 * th_11
    z1 = N1_10 * th_10 + N1_11 * th_11
    x2 = N2_00 * th_20 + N2_01 * th_21
    z2 = N2_10 * th_20 + N2_11 * th_21

    Q1 = E1([x1, z1])
    Q2 = E2([x2, z2])

    return CoupleECPoint(Q1, Q2)


def splitting_eval(pts, E1, E2, N1, N2, hadamard):
    r"""
    [Algorithm 4.130] SplittingEval.

    INPUT:

    - ``pts`` -- a list of :class:`DimTwoPoint`
    - ``E1``, ``E2`` -- the curves of the product
    - ``N1``, ``N2`` -- the conversion matrices for ``E1`` and ``E2``
    - ``hadamard`` -- a boolean

    OUTPUT:

    The list of the points of ``pts`` as :class:`CoupleECPoint`.
    """
    image_pts = []
    for R in pts:
        P = dim_two_point_to_ec_points(R, E1, E2, N1, N2, hadamard)
        image_pts.append(P)
    return image_pts
