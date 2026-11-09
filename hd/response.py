from sage.all import gcd, inverse_mod

import quat.ideals as qt_ids
import id2iso.id2iso as qlpt
from hd.isogeny_chain import isogeny_two_two_chain
from hd.couple_point import CoupleECPoint
from ec.ec_basis import ECBasis
import crypto.utils as utils

def hd_response(I_sk_chl, alpha_rsp, alpha_sk_chl, q_rsp, E_com, B_com,
                prng_res, prng_def, ec_ctx, quat_ctx):
    r"""
    [Algorithm 3.21] HDResponse.

    Compute the auxiliary and challenge curves of the response, with their
    torsion bases.

    INPUT:

    - ``I_sk_chl`` -- an ideal
    - ``alpha_rsp`` -- an element of ``I_sk_chl`` of norm less than
      `2^{e_{rsp} - 1}`
    - ``alpha_sk_chl`` -- an element of ``I_sk_chl``
    - ``q_rsp`` -- an integer
    - ``E_com`` -- the commitment curve
    - ``B_com`` -- a torsion basis of ``E_com``
    - ``prng_res`` -- the response pseudorandom generator
    - ``prng_def`` -- the default-domain pseudorandom generator
    - ``ec_ctx`` -- the elliptic curve context
    - ``quat_ctx`` -- the quaternion context

    OUTPUT:

    A tuple ``(E_aux, B_aux, E_chl, B_chl)`` with the auxiliary and challenge
    curves and their `2^{e_{rsp} + 2}`-torsion bases.
    """
    e_rsp = ec_ctx.e_rsp

    if gcd(2**e_rsp - q_rsp, alpha_sk_chl.reduced_norm()) != 1:
        raise ValueError("Auxiliary degree not comprime to sk/chall degree")

    I_aux, theta_aux = qt_ids.random_ideal_given_norm(
            2**e_rsp - q_rsp, alpha_sk_chl, prng_res, quat_ctx)


    # Inert intersection
    I_asc = ((I_aux * theta_aux).intersection(I_sk_chl)) * theta_aux**-1

    E_aux, B_aux = qlpt.ideal_to_isogeny(I_asc, prng_def, quat_ctx, ec_ctx)

    f = quat_ctx.f
    M_rsp = utils.endomorphism_to_matrix(ec_ctx, alpha_rsp)
    M_aux = utils.endomorphism_to_matrix(ec_ctx, theta_aux)
    mu = inverse_mod(I_sk_chl.norm(), 2**f)
    B_aux = B_aux.act_on_basis(mu * M_aux * M_rsp)

    q_inv = inverse_mod(q_rsp, 2**(e_rsp + 2))

    P_com, Q_com, PQ_com = [P.ec_ladder(2**(f - e_rsp - 2)) for P in B_com]
    P_aux, Q_aux, _ = [P.ec_ladder(q_inv * 2**(f - e_rsp - 2)) for P in B_aux]

    T1 = CoupleECPoint(P_com, P_aux)
    T2 = CoupleECPoint(Q_com, Q_aux)
    X1 = CoupleECPoint(P_com, E_aux.zero())
    X2 = CoupleECPoint(Q_com, E_aux.zero())
    X12 = CoupleECPoint(PQ_com, E_aux.zero())
    E_aux, E_chl, [Y1, Y2, Y12] = isogeny_two_two_chain(
            E_com, E_aux, T1, T2, e_rsp, [X1, X2, X12], False)

    P_aux = Y1.P1
    Q_aux = Y2.P1
    PQ_aux = Y12.P1
    B_aux = ECBasis(P_aux, Q_aux, PQ_aux)

    P_chl = Y1.P2
    Q_chl = Y2.P2
    PQ_chl = Y12.P2
    B_chl = ECBasis(P_chl, Q_chl, PQ_chl)

    return E_aux, B_aux, E_chl, B_chl

