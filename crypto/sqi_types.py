import crypto.encode as benc
from ec.ec_curve import ECCurve

class PublicKey:
    r"""
    Implements the type PublicKey [Section 3.1.3].

    INPUT:

    - ``Epk`` -- the public curve, an :class:`ECCurve`
    - ``hint`` -- the hint ``(h, hA)`` of the torsion basis on ``Epk``
    """

    def __init__(self, Epk, hint):
        self.Epk = Epk
        self.hint = hint

    @classmethod
    def from_hex(cls, pk_hex, ec_ctx):
        r"""
        Decode a public key with :func:`hex_to_pk`.

        INPUT:

        - ``pk_hex`` -- a hex string encoding ``PUBLICKEY_BYTES`` bytes
        - ``ec_ctx`` -- the elliptic curve context
        """
        A, hint = benc.hex_to_pk(pk_hex, ec_ctx)
        return cls.from_ints(A, hint, ec_ctx)

    @classmethod
    def from_ints(cls, A, hint, ec_ctx):
        r"""
        Build a public key from integer data.

        INPUT:

        - ``A`` -- the integers ``(a, b)`` of the coefficient `a + bi` of
          ``Epk``
        - ``hint`` -- the hint ``(h, hA)``
        - ``ec_ctx`` -- the elliptic curve context
        """
        Epk = ECCurve(ec_ctx.Fp2(list(A)))
        return cls(Epk, hint)

    def __eq__(self, other):
        if not isinstance(other, PublicKey):
            return False
        return self.Epk == other.Epk and self.hint == other.hint

    def __repr__(self):
        return f"Public Key with curve {self.Epk}"

class SecretKey:
    r"""
    Implements the type SecretKey [Section 3.1.3].

    INPUT:

    - ``pk`` -- the associated :class:`PublicKey`
    - ``Isk`` -- the inert form ``(N, x, y)`` of the secret ideal
    - ``M`` -- the matrix ``(a, b, c, d)`` mapping the deterministic basis on
      `E_{pk}` to the secret basis
    """

    def __init__(self, pk, Isk, M):
        self.pk = pk
        self.Isk = Isk
        self.M = M

    @classmethod
    def from_hex(cls, sk_hex, ec_ctx):
        r"""
        Decode a secret key with :func:`hex_to_sk`.

        INPUT:

        - ``sk_hex`` -- a hex string encoding ``SECRETKEY_BYTES`` bytes
        - ``ec_ctx`` -- the elliptic curve context
        """
        pk, Isk, M = benc.hex_to_sk(sk_hex, ec_ctx)
        pk = PublicKey.from_ints(*pk, ec_ctx)
        return cls(pk, Isk, M)

    def __eq__(self, other):
        if not isinstance(other, SecretKey):
            return False
        return self.pk == other.pk and self.Isk == other.Isk and self.M == other.M

    def __repr__(self):
        N, x, y = self.Isk
        return f"Secret Key with ideal (N, x, y) = ({N}, {x}, {y})"

class Signature:
    r"""
    Implements the type Signature [Section 3.1.3], together with the signed
    message.

    INPUT:

    - ``chl`` -- the challenge, an integer
    - ``hint_chl`` -- the hint ``(h, hA)`` of the basis on the challenge curve
    - ``M_chl`` -- the matrix ``(a, b, c, d)`` acting on the deterministic
      basis of the challenge curve
    - ``E_aux`` -- the auxiliary curve, an :class:`ECCurve`
    - ``hint_aux`` -- the hint ``(h, hA)`` of the basis on ``E_aux``
    - ``msg`` -- the message, as bytes
    - ``ec_ctx`` -- the elliptic curve context
    """

    def __init__(self, chl, hint_chl, M_chl, E_aux, hint_aux, msg, ec_ctx):
        self.chl = chl
        self.hint_chl = hint_chl
        self.M_chl = M_chl
        self.E_aux = E_aux
        self.hint_aux = hint_aux
        self.msg = msg
        self.ec_ctx = ec_ctx

    @classmethod
    def from_hex(cls, sig_hex, mlen, ec_ctx):
        r"""
        Decode a signed message with :func:`hex_to_sig`.

        INPUT:

        - ``sig_hex`` -- a hex string encoding ``SIGNATURE_BYTES + mlen``
          bytes
        - ``mlen`` -- the message length in bytes
        - ``ec_ctx`` -- the elliptic curve context
        """
        chl, hint_chl, M_chl, A_aux, hint_aux, msg = benc.hex_to_sig(
                sig_hex, mlen, ec_ctx)
        E_aux = ECCurve(ec_ctx.Fp2(list(A_aux)))
        return cls(chl, hint_chl, M_chl, E_aux, hint_aux, msg, ec_ctx)

    def __eq__(self, other):
        if not isinstance(other, Signature):
            return False

        return all([
            self.chl == other.chl,
            self.hint_chl == other.hint_chl,
            self.M_chl == other.M_chl,
            self.E_aux == other.E_aux,
            self.hint_aux == other.hint_aux,
            self.msg == other.msg,
        ])

    def __repr__(self):
        return f"Signature for message {self.msg.hex()}"

