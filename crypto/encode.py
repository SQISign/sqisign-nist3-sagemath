import hashlib

# Hash
def hash_to_challenge(pk, Ecom, msg, ec_ctx):
    r"""
    [Algorithm 3.3] HASH.

    INPUT:

    - ``pk`` -- the public key
    - ``Ecom`` -- the commitment curve
    - ``msg`` -- the message, as bytes of arbitrary length
    - ``ec_ctx`` -- the elliptic curve context

    OUTPUT:

    The challenge, as the little-endian integer of the ``CHALLENGE_BYTES``
    bytes of output.
    """
    shake = hashlib.shake_256()
    shake.update(b'SQI')
    buf = pk_to_bytes(pk, ec_ctx) + curve_to_bytes(Ecom, ec_ctx)
    shake.update(buf)
    shake.update(msg)
    out = shake.digest(ec_ctx.lb // 8)
    return int.from_bytes(out, byteorder='little')

# Encoders
def fp2_to_bytes(x, ec_ctx):
    r"""
    Encode an element of `\mathbb{F}_{p^2}` [Chapter 6].

    INPUT:

    - ``x`` -- an element of `\mathbb{F}_{p^2}`
    - ``ec_ctx`` -- the elliptic curve context

    OUTPUT:

    ``2*FP_ENCODED_BYTES`` bytes: the real part followed by the imaginary part,
    each as a little-endian integer of ``FP_ENCODED_BYTES`` bytes.
    """
    a, b = x
    a_bytes = int.to_bytes(int(a), ec_ctx.FP_ENCODED_BYTES, byteorder='little')
    b_bytes = int.to_bytes(int(b), ec_ctx.FP_ENCODED_BYTES, byteorder='little')

    return a_bytes + b_bytes

def curve_to_bytes(x, ec_ctx):
    r"""
    Encode a curve by its Montgomery coefficient [Chapter 6].

    INPUT:

    - ``x`` -- an :class:`ECCurve`
    - ``ec_ctx`` -- the elliptic curve context

    OUTPUT:

    ``2*FP_ENCODED_BYTES`` bytes.
    """
    return fp2_to_bytes(x.A_coefficient(), ec_ctx)

def pk_to_bytes(pk, ec_ctx):
    r"""
    Encode a public key [Chapter 6].

    INPUT:

    - ``pk`` -- the public key
    - ``ec_ctx`` -- the elliptic curve context

    OUTPUT:

    ``PUBLICKEY_BYTES = 2*FP_ENCODED_BYTES + 1`` bytes: the curve `E_{pk}`
    followed by the hint byte `h_A + 2h`.
    """
    h, ha = pk.hint
    return curve_to_bytes(pk.Epk, ec_ctx) + bytes([h * 2 + ha])

# Decoders
def hex_to_fp2(x, ec_ctx):
    r"""
    Decode an element of `\mathbb{F}_{p^2}` [Chapter 6].

    INPUT:

    - ``x`` -- a hex string encoding ``2*FP_ENCODED_BYTES`` bytes
    - ``ec_ctx`` -- the elliptic curve context

    OUTPUT:

    A tuple ``(a, b)`` of integers representing `a + bi`.
    """
    ab = bytes.fromhex(x)
    if not len(ab) == 2 * ec_ctx.FP_ENCODED_BYTES:
        raise ValueError("Malformed input to hex_to_fp2")

    a = int.from_bytes(ab[:ec_ctx.FP_ENCODED_BYTES], byteorder='little')
    b = int.from_bytes(ab[ec_ctx.FP_ENCODED_BYTES:], byteorder='little')

    return (a, b)

def hex_to_hint(x):
    r"""
    Decode a hint byte [Chapter 6].

    INPUT:

    - ``x`` -- an integer in `[0, 256)`, the hint byte `h_A + 2h`

    OUTPUT:

    The hint ``(h, hA)``.
    """
    return (x // 2, x % 2)

def hex_to_pk(x, ec_ctx):
    r"""
    Decode a public key [Chapter 6].

    INPUT:

    - ``x`` -- a hex string encoding ``PUBLICKEY_BYTES = 2*FP_ENCODED_BYTES + 1``
      bytes
    - ``ec_ctx`` -- the elliptic curve context

    OUTPUT:

    A tuple ``(A, hint)`` with ``A = (a, b)`` the Montgomery coefficient of
    `E_{pk}` and ``hint = (h, hA)``.
    """
    A_len = 4 * ec_ctx.FP_ENCODED_BYTES
    if not len(x) == A_len + 2:
        raise ValueError("Malformed input to hex_to_pk")
    A = hex_to_fp2(x[:A_len], ec_ctx)

    hint = hex_to_hint(int.from_bytes(bytes.fromhex(x[-2:]), byteorder='little'))
    return (A, hint)

def hex_to_ideal(x, ec_ctx):
    r"""
    Decode an inert ideal [Chapter 6].

    INPUT:

    - ``x`` -- a hex string encoding ``3*FP_ENCODED_BYTES`` bytes
    - ``ec_ctx`` -- the elliptic curve context

    OUTPUT:

    The inert form ``(N, x, y)`` of the ideal.
    """
    FPB = ec_ctx.FP_ENCODED_BYTES
    if not len(x) == 6 * FPB:
        raise ValueError("Malformed input to hex_to_ideal")
    x = bytes.fromhex(x)

    N = int.from_bytes(x[:FPB], byteorder='little')
    x_id = int.from_bytes(x[FPB:2*FPB], byteorder='little')
    y_id = int.from_bytes(x[2*FPB:], byteorder='little')
    return (N, x_id, y_id)

def hex_to_matrix(x, n_bytes):
    r"""
    Decode a `2 \times 2` integer matrix [Chapter 6].

    INPUT:

    - ``x`` -- a hex string encoding ``4*n_bytes`` bytes
    - ``n_bytes`` -- the number of bytes of each entry

    OUTPUT:

    The entries ``(a, b, c, d)`` of the matrix, in row order.
    """
    if not len(x) == 8 * n_bytes:
        raise ValueError("Malformed input to hex_to_matrix")
    x = bytes.fromhex(x)

    return tuple(int.from_bytes(x[i*n_bytes:(i+1)*n_bytes], byteorder='little')
                 for i in range(4))

def hex_to_sk(x, ec_ctx):
    r"""
    Decode a secret key [Chapter 6].

    INPUT:

    - ``x`` -- a hex string encoding ``SECRETKEY_BYTES = 5*FP_ENCODED_BYTES +
      4*CHALLENGE_BYTES + 1`` bytes: the public key
      (``2*FP_ENCODED_BYTES + 1``), the secret ideal (``3*FP_ENCODED_BYTES``)
      and the matrix (``4*CHALLENGE_BYTES``)
    - ``ec_ctx`` -- the elliptic curve context

    OUTPUT:

    A tuple ``(pk, Isk, M)`` with ``pk`` as in :func:`hex_to_pk`,
    ``Isk = (N, x, y)`` and ``M = (a, b, c, d)``.
    """
    pk_len = 4 * ec_ctx.FP_ENCODED_BYTES + 2
    Isk_len = 6 * ec_ctx.FP_ENCODED_BYTES
    M_len = 8 * ec_ctx.CHALLENGE_BYTES
    if not len(x) == pk_len + Isk_len + M_len:
        raise ValueError("Malformed input to hex_to_sk")

    pk, x = x[:pk_len], x[pk_len:]
    Isk, M = x[:Isk_len], x[Isk_len:]

    pk = hex_to_pk(pk, ec_ctx)
    Isk = hex_to_ideal(Isk, ec_ctx)
    M = hex_to_matrix(M, ec_ctx.CHALLENGE_BYTES)
    return (pk, Isk, M)

def hex_to_sig(x, mlen, ec_ctx):
    r"""
    Decode a signed message, i.e. a signature followed by the message
    [Chapter 6].

    INPUT:

    - ``x`` -- a hex string encoding ``SIGNATURE_BYTES + mlen`` bytes, with
      ``SIGNATURE_BYTES = 2*FP_ENCODED_BYTES + 4*RESPONSE_BYTES +
      CHALLENGE_BYTES + 2``: the curve `E_{aux}` (``2*FP_ENCODED_BYTES``), the
      matrix `M_{chl}` (``4*RESPONSE_BYTES``), the challenge
      (``CHALLENGE_BYTES``), the hints of `E_{aux}` and `E_{chl}` (one byte
      each) and the message (``mlen``)
    - ``mlen`` -- the message length in bytes
    - ``ec_ctx`` -- the elliptic curve context

    OUTPUT:

    A tuple ``(chl, hint_chl, M_chl, A_aux, hint_aux, msg)``.
    """
    RESP_BYTES = ec_ctx.RESPONSE_BYTES
    CHALL_BYTES = ec_ctx.CHALLENGE_BYTES
    FP_BYTES = ec_ctx.FP_ENCODED_BYTES

    if not len(x) == 2 * (2*FP_BYTES + 4*RESP_BYTES + CHALL_BYTES + 2 + mlen):
        raise ValueError("Malformed input to hex_to_sig")

    A_aux, x = x[:4*FP_BYTES], x[4*FP_BYTES:]
    A_aux = hex_to_fp2(A_aux, ec_ctx)

    M_chl, x = x[:8*RESP_BYTES], x[8*RESP_BYTES:]
    M_chl = hex_to_matrix(M_chl, RESP_BYTES)

    x = bytes.fromhex(x)

    chl = int.from_bytes(x[:CHALL_BYTES], byteorder='little')
    x = x[CHALL_BYTES:]

    hint_aux = hex_to_hint(x[0])
    hint_chl = hex_to_hint(x[1])
    msg = x[2:]

    return (chl, hint_chl, M_chl, A_aux, hint_aux, msg)
