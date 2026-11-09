import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import sqisign as sqi
from crypto.aes256_ctr_drbg import AES256_CTR_DRBG

MSG_BYTES = 32


def benchmark(lvl=1, n_runs=10, seed=None):
    """
    Time keygen, signing and verification n_runs times on random messages
    """
    if seed is None:
        seed = os.urandom(48).hex()
    print(f"Level {lvl}, {n_runs} runs, seed {seed}")

    # Keys, signatures and messages are all derived from the seed
    drbg = AES256_CTR_DRBG(bytes.fromhex(seed))
    SQI = sqi.SQIsign(level=lvl, randombytes=drbg.random_bytes)

    t_keygen, t_sign, t_verify = 0, 0, 0
    for k in range(n_runs):
        msg = drbg.random_bytes(MSG_BYTES)

        t0 = time.time()
        SQI.keygen()
        t1 = time.time()
        sig = SQI.sign(msg)
        t2 = time.time()
        valid = sqi.verify(SQI.pk, sig)
        t3 = time.time()

        if not valid:
            raise RuntimeError(f"[Run {k}] signature does not verify")

        t_keygen += t1 - t0
        t_sign += t2 - t1
        t_verify += t3 - t2
        print(f"Run {k}: keygen {t1 - t0:.2f}s, "
              f"sign {t2 - t1:.2f}s, verify {t3 - t2:.2f}s")

    print(f"Average: keygen {t_keygen / n_runs:.2f}s, "
          f"sign {t_sign / n_runs:.2f}s, verify {t_verify / n_runs:.2f}s")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Benchmark SQIsign keygen, signing and verification."
    )
    parser.add_argument("-l", "--level", type=int, choices=[1, 3, 5],
                        default=1, help="security level (default: 1)")
    parser.add_argument("-n", "--n_runs", type=int, default=10,
                        help="number of runs (default: 10)")
    parser.add_argument("-s", "--seed", type=str, default=None,
                        help="48-byte DRBG seed as hex (default: random)")
    args = parser.parse_args()

    benchmark(args.level, args.n_runs, args.seed)
