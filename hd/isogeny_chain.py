from hd.gluing_isogeny import first_isogeny
from hd.theta import (
    compute_constants,
    compute_constants_from_iso_data,
    recover_dual_null,
    generic_codomain,
    generic_eval,
    dim_two_dbl_iter,
)
from hd.splitting import dim_two_struct_to_ec_curves, splitting_eval


def isogeny_two_two_chain(E1, E2, T1, T2, e, pts, verify):
    r"""
    [Algorithm 4.132] Isogeny22Chain.

    Compute the `(2, 2)`-isogeny chain `\Phi : E_1 \times E_2 \to E_3 \times
    E_4` of kernel `\langle [4]T_1, [4]T_2 \rangle` and evaluate it on ``pts``.

    INPUT:

    - ``E1``, ``E2`` -- the domain curves
    - ``T1``, ``T2`` -- points of type :class:`CoupleECPoint`
    - ``e`` -- the length of the chain
    - ``pts`` -- a list of :class:`CoupleECPoint` of the form `(P_1, 0)` or
      `(0, P_2)`
    - ``verify`` -- if ``True``, verify the kernel at every step

    OUTPUT:

    A tuple ``(E3, E4, image_pts)`` with the codomain curves and the images
    of ``pts``.
    """
    # Bookkeeping for values used when m != 1 in the loop
    x1, x2, y1, y2 = (None, None, None, None)

    B, orders, theta_q1, theta_q2, image_points, k = first_isogeny(
        E1, E2, T1, T2, e, pts, verify
    )
    U, V = compute_constants(B)

    for m in range(1, e):
        not_step_two = m != 1
        if not_step_two:
            B, V = compute_constants_from_iso_data(U, x1, x2, y1, y2)

        while orders[k] != 1:
            assert orders[k] >= 2
            assert (k + 1) < len(orders)

            n = orders[k] // 2
            theta_q1[k + 1] = dim_two_dbl_iter(theta_q1[k], U, V, n)
            theta_q2[k + 1] = dim_two_dbl_iter(theta_q2[k], U, V, n)
            orders[k + 1] = orders[k] - n
            k = k + 1

        K1 = theta_q1[k]
        K2 = theta_q2[k]

        # Compute the (2,2)-isogeny, for the last step we delay verification
        if m == e - 1:
            U, x1, x2, y1, y2 = generic_codomain(K1, K2, True, False)
        else:
            U, x1, x2, y1, y2 = generic_codomain(K1, K2, not_step_two, verify)

        # Push image points through
        image_points = [generic_eval(P, U, not_step_two) for P in image_points]

        # Push the kernel through
        assert orders[k] == 1
        for i in range(k):
            theta_q1[i] = generic_eval(theta_q1[i], U, not_step_two)
            theta_q2[i] = generic_eval(theta_q2[i], U, not_step_two)
            orders[i] -= 1
        k -= 1

    # Recover B from U
    B = recover_dual_null(U, x1, x2, y1, y2)

    # Recover the elliptic curves and isomorphisms from B
    E3, E4, N3, N4 = dim_two_struct_to_ec_curves(B, True)

    # Recover the image points
    image_points = splitting_eval(image_points, E3, E4, N3, N4, True)

    return E3, E4, image_points
