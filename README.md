# SQIsign SageMath

A SageMath implementation of SQIsign following the SQIsign specification
(round 3).

## Requirements

- [SageMath](https://www.sagemath.org/)
- [pycryptodome](https://pypi.org/project/pycryptodome/), install with:

  ```
  sage -pip install pycryptodome
  ```

## Project overview

- `sqisign.py`: key generation, signing and verification
- `quat/`: quaternions and ideals
- `id2iso/`: ideal-to-isogeny translation
- `ec/`: Montgomery curves and isogenies
- `hd/`: (2,2)-isogenies in the theta model
- `crypto/`: PRNG, encodings, bounded integers and utilities
- `params/`: precomputed data for each security level
- `KAT/`: known answer tests from the C implementation
- `tests/`: KAT validation and benchmarks

## Usage

To sign a message and verify a signature:

```python
from sqisign import SQIsign, verify

sqi = SQIsign(level=1)
sqi.keygen()
sig = sqi.sign(b"Hello SQIsign")
assert verify(sqi.pk, sig)
```

### Benchmarking

To benchmark the code run:

```
sage tests/benchmark.py
```

which will average over 10 runs of keygen, signing and verification at level 1.
Additional options are:
- `-l level`, 1/3/5, default 1
- `-n n_runs`, number of runs, default 10
- `-s seed`, a 48-byte DRBG seed (default random)

### KAT verification

To test that the code passes all the KATs, run

```
sage tests/kat_verify.py
```

Additional options are:
- `-l level`, 1/3/5, default 1
- `-k kat_number`, to validate a specific KAT instance; by default all the
  instances of the chosen level are validated
If no options are provided the code tests all the KATs.

## A note on randomness

All randomness is drawn from the domain-separated PRNGs of the
specification (`"key"`, `"com"`, `"res"` and `"def"`), seeded from
`randombytes`. By default `randombytes` is `os.urandom`; the KAT and benchmark
scripts pass a seeded `AES256_CTR_DRBG` instead.

According to the specification, the `"def"` PRNG only feeds Las Vegas
subroutines, whose output does not depend on the randomness, so it can be
replaced by any other source. This is not the case at level 5: there the
prime `Dmix` is `1 mod 4`, and the square root computed in
`prime_norm_inert_ideal_sampling` uses Tonelli-Shanks, whose output sign
depends on the non-residue sampled from `"def"`. Reproducing the level 5 KATs
therefore requires the `"def"` stream to match the C implementation exactly. In
all other cases the `"def"` sream can be replaced by any other stream
(including `sage` builtin prng, as already done in some algorithms) without
altering the result.
