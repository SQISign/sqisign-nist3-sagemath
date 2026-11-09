from Crypto.Hash import SHAKE256

import crypto.bint as bint

PRNG_SEED_BYTES = 48
PRNG_DOMAIN_SEP_BYTES = 3
DS_DEFAULT_DOMAIN = b"def"

class PRNG:
    r"""
    Implements the type PRNG [Section 3.3.1]

    The stream SHAKE256(``root_seed`` || ``domain``), squeezed incrementally.

    INPUT:

    - ``root_seed`` -- a 48-byte seed
    - ``domain`` -- a 3-byte domain separator, as bytes or an ASCII string
      (default: ``b"def"``)
    """

    def __init__(self, root_seed, domain=DS_DEFAULT_DOMAIN):
        root_seed = bytes(root_seed)
        if len(root_seed) != PRNG_SEED_BYTES:
            raise ValueError(
                f"root_seed must be {PRNG_SEED_BYTES} bytes, "
                f"got {len(root_seed)}"
            )
        if isinstance(domain, str):
            domain = domain.encode("ascii")
        if len(domain) != PRNG_DOMAIN_SEP_BYTES:
            raise ValueError(
                f"domain separator must be {PRNG_DOMAIN_SEP_BYTES} bytes, "
                f"got {len(domain)}"
            )
        self.domain = domain
        self._xof = SHAKE256.new(root_seed + domain)

    def random_bytes(self, n):
        r"""
        [Algorithm 3.2] PRNG.Bytes.

        INPUT:

        - ``n`` -- a non-negative integer

        OUTPUT:

        The next ``n`` pseudorandom bytes of the stream.
        """
        return self._xof.read(n)

    def sample_from_interval(self, a, b):
        r"""
        [Algorithm 4.24] SampleFromInterval.

        Sample a pseudorandom integer in `[a, b]`, using `\mathrm{bound}(b) + 64`
        random bits.

        INPUT:

        - ``a`` -- an integer or :class:`BoundedInt` with `-b \le a < b`
        - ``b`` -- a non-negative integer or :class:`BoundedInt`

        OUTPUT:

        An integer in `[a, b]`.
        """
        if not isinstance(b, bint.BoundedInt):
            b = bint.BoundedInt(b)

        if isinstance(a, bint.BoundedInt):
            a = a.nn

        m = b.nn - a + 1
        n = b.bb + 64

        nbytes = (n + 7) // 8
        mask = (1 << n) - 1

        while True:
            buf = self.random_bytes(nbytes)
            x = (int.from_bytes(buf, "little") & mask) + m

            if x < (1 << n):
                break

        return a + (x % m)

