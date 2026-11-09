import os

from sage.all import ZZ, Zmod, matrix

from quat.quat_ctx import QuatContext
from ec.ec_ctx import ECContext
from crypto.prng import PRNG, PRNG_SEED_BYTES

from crypto.sqi_types import PublicKey, SecretKey, Signature

import crypto.utils as utils
import id2iso.id2iso as qlpt
import quat.ideals as qt_ids
import quat.inert_ideals as qt_inert
import crypto.encode as benc
from ec.ec_isogeny import ECIsogeny
from hd.couple_point import CoupleECPoint
from hd.isogeny_chain import isogeny_two_two_chain
import hd.response as hdresp

class SQIsign:
    r"""
    SQIsign signer at a given NIST security level, holding a key pair.

    INPUT:

    - ``level`` -- the NIST level, one of ``1``, ``3`` or ``5`` (default: ``1``)
    - ``randombytes`` -- a callable ``n -> n`` random bytes, the analogue of the
      C ``randombytes()`` (default: ``os.urandom``); pass a seeded
      ``AES256_CTR_DRBG.random_bytes`` for KAT reproduction
    """

    def __init__(self, level=1, randombytes=None):
        self.randombytes = randombytes if randombytes is not None \
            else os.urandom

        self.sk = None
        self.pk = None

        if level == 1:
            self.p = 3 * 2**324 - 1
            self.lb = 128
            self.quat_ctx = QuatContext(self.p, self.lb)
            self.ec_ctx = ECContext(self.p, self.lb)

        elif level == 3:
            self.p = 27 * 2**500 - 1
            self.lb = 192
            self.quat_ctx = QuatContext(self.p, self.lb)
            self.ec_ctx = ECContext(self.p, self.lb)

        elif level == 5:
            self.p = 17 * 2**664 - 1
            self.lb = 256
            self.quat_ctx = QuatContext(self.p, self.lb)
            self.ec_ctx = ECContext(self.p, self.lb)

        else:
            raise ValueError(f"level {level} not recognised")

    def seed_prng(self):
        r"""
        Draw a fresh 48-byte seed from ``randombytes`` and initialize the
        ``"def"``, ``"key"``, ``"com"`` and ``"res"`` PRNGs from it, as in the
        first steps of [Algorithm 3.22] and [Algorithm 3.23].
        """
        self.root_seed = self.randombytes(PRNG_SEED_BYTES)

        # Default domain (only relevant at level 5)
        self.prng_def = PRNG(self.root_seed, "def")

        # Key generation randomness
        self.prng_key = PRNG(self.root_seed, "key")

        # Signing randomness (commitment - response)
        self.prng_com = PRNG(self.root_seed, "com")
        self.prng_res = PRNG(self.root_seed, "res")

    def keygen(self):
        r"""
        [Algorithm 3.22] SQIsign.KeyGen.

        """
        self.seed_prng()

        I_sk = qt_ids.uniform_O0_class_sampling(0, self.prng_key, self.prng_def, self.quat_ctx)
        E_pk, phi_sk = qlpt.ideal_to_isogeny(
            I_sk, self.prng_def, self.quat_ctx, self.ec_ctx)

        Bpk, hint_pk = E_pk.torsion_basis_to_hint(self.lb+2, 0)
        Msk = phi_sk.change_of_basis(Bpk, self.lb+2)
        Msk = tuple(x % 2**self.lb for x in Msk.list())
        Msk = utils.normalize_matrix(Msk, self.lb)

        I_sk = qt_inert.ideal_to_Nxy(I_sk)

        self.pk = PublicKey(E_pk, hint_pk)
        self.sk = SecretKey(self.pk, I_sk, Msk)

    def sign(self, msg):
        r"""
        [Algorithm 3.23] SQIsign.Sign.

        INPUT:

        - ``msg`` -- the message, as bytes

        OUTPUT:

        A :class:`Signature`.
        """
        self.seed_prng() # Reinstantiate randomness

        if self.pk == None or self.sk == None:
            self.keygen()

        sk = self.sk

        # Commitment
        I_com = qt_ids.uniform_O0_class_sampling(2, self.prng_com, self.prng_def, self.quat_ctx)
        E_com, B_com = qlpt.ideal_to_isogeny(
            I_com, self.prng_def, self.quat_ctx, self.ec_ctx)

        # Challenge
        chl = benc.hash_to_challenge(sk.pk, E_com, msg, self.ec_ctx)

        # Response
        c1 = (sk.M[0] + sk.M[1] * chl) % 2**self.lb
        c2 = (sk.M[2] + sk.M[3] * chl) % 2**self.lb

        I_chl, gamma = qt_ids.kernel_to_ideal(
                self.lb, c1, c2, self.ec_ctx, self.quat_ctx)
        I_sk = qt_inert.xy_to_ideal(*sk.Isk, self.quat_ctx)

        # Inert intersection
        I_sk_chl = ((I_chl * gamma).intersection(I_sk)) * gamma**-1

        I_sk_chl, alpha_rsp, alpha_sk_chl, q_rsp = qt_ids.sample_quat_response(
                I_sk_chl, I_com, self.prng_res, self.prng_def, self.quat_ctx,
                self.ec_ctx)

        # Response Packaging
        E_aux, B_aux, E_chl, B_chl = hdresp.hd_response(
            I_sk_chl, alpha_rsp, alpha_sk_chl, q_rsp, E_com, B_com,
            self.prng_res, self.prng_def, self.ec_ctx, self.quat_ctx)

        M_chl, hint_aux, hint_chl = utils.set_change_of_basis_matrix(
                E_aux, E_chl, B_aux, B_chl, self.ec_ctx.e_rsp + 2)
        M_chl = utils.normalize_matrix(
                tuple(ZZ(x) for x in M_chl.list()), self.ec_ctx.e_rsp + 2)

        sig = Signature(chl, hint_chl, M_chl, E_aux, hint_aux, msg, self.ec_ctx)
        return sig


def verify(pk, sig):
    r"""
    [Algorithm 3.24] SQIsign.Verify.

    INPUT:

    - ``pk`` -- a :class:`PublicKey`
    - ``sig`` -- a :class:`Signature`, carrying the signed message

    OUTPUT:

    ``True`` if ``sig`` is a valid signature of ``sig.msg`` under ``pk``,
    ``False`` otherwise.
    """
    # Loading context
    ec_ctx = sig.ec_ctx
    e_rsp = ec_ctx.e_rsp
    Epk = pk.Epk

    # Sanity checks
    if not all(0 <= x < 2**(e_rsp+2) for x in sig.M_chl):
        return False

    if not sig.M_chl == utils.normalize_matrix(sig.M_chl, e_rsp+2):
        return False

    # Challenge kernel on Epk
    Ppk, Qpk, PQpk = Epk.torsion_basis_from_hint(ec_ctx.lb + 2, pk.hint)
    Kchl = Ppk.ec_ladder3pt(Qpk, PQpk, sig.chl)

    # Compute Echl from the kernel
    phi_chl = ECIsogeny(Epk, Kchl, ec_ctx.lb)
    Echl, _ = phi_chl.ec_two_isogeny_chain([])

    # Deterministic basis on Eaux of order e_rsp+2
    Eaux = sig.E_aux
    Paux, Qaux, _ = Eaux.torsion_basis_from_hint(e_rsp + 2, sig.hint_aux)

    # Challenge basis on Echl
    Echl_bs = Echl.torsion_basis_from_hint(e_rsp + 2, sig.hint_chl)
    M_chl = matrix(Zmod(2**(e_rsp + 2)), 2, sig.M_chl)
    Pchl, Qchl, _ = Echl_bs.act_on_basis(M_chl)

    # Compute codomain with 2D isogenies
    T1 = CoupleECPoint(Pchl, Paux)
    T2 = CoupleECPoint(Qchl, Qaux)

    try:
        _, Ecom, _ = isogeny_two_two_chain(Echl, Eaux, T1, T2, e_rsp, [], False)
    except ValueError:
        return False

    # Recover the challenge from the commitment curve
    rec_chl = benc.hash_to_challenge(pk, Ecom, sig.msg, ec_ctx)

    return rec_chl == sig.chl
