import sys
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
KAT_DIR = ROOT / "KAT"
sys.path.insert(0, str(ROOT))

import sqisign as sqi
from crypto.aes256_ctr_drbg import AES256_CTR_DRBG

from crypto.sqi_types import PublicKey, SecretKey, Signature

def validate(lvl=1, k_test=None):
    """
    Validate KATs for the given level
    """
    if lvl == 1:
        filename = KAT_DIR / "PQCsignKAT_270_SQIsign_p324_3"
    elif lvl == 3:
        filename = KAT_DIR / "PQCsignKAT_417_SQIsign_p500_27"
    elif lvl == 5:
        filename = KAT_DIR / "PQCsignKAT_549_SQIsign_p664_17"
    else:
        raise ValueError(f"Invalid level: {lvl}")

    # Patterns
    kat_pattern = re.compile(
        r"count = (\d+)\n"
        r"seed = ([0-9A-Fa-f]+)\n"
        r"mlen = (\d+)\n"
        r"msg = ([0-9A-Fa-f]+)\n"
        r"pk =\s*([0-9A-Fa-f]*)\n"
        r"sk =\s*([0-9A-Fa-f]*)\n"
        r"smlen =\s*(\d*)\n"
        r"sm =\s*([0-9A-Fa-f]*)"
    )

    req = {}
    rsp = {}

    # Parse req
    with open(f"{filename}.req", "r") as f:
        req_content = f.read()

    blocks = [b for b in req_content.split("\n\n") if b.strip() and not b.startswith("#")]

    for block in blocks:
        match = kat_pattern.search(block)
        assert match is not None
        count = int(match.group(1))
        data = {
            "seed": match.group(2),
            "mlen": int(match.group(3)),
            "msg": match.group(4),
        }
        req[count] = data

    # Parse rsp
    with open(f"{filename}.rsp", "r") as f:
        rsp_content = f.read()

    rsp_blocks = [b for b in rsp_content.split("\n\n") if b.strip() and not b.startswith("#")]

    for block in rsp_blocks:
        match = kat_pattern.search(block)
        assert match is not None
        count = int(match.group(1))
        data = {
            "seed": match.group(2),
            "mlen": int(match.group(3)),
            "msg": match.group(4),
            "pk": match.group(5),
            "sk": match.group(6),
            "smlen": int(match.group(7)) if match.group(7) else "",
            "sm": match.group(8),
        }
        rsp[count] = data

    if k_test is not None and k_test not in req:
        raise ValueError(f"Invalid KAT number for level {lvl}: {k_test} (available: {min(req)}-{max(req)})")

    # Validate KAT
    SQI = sqi.SQIsign(level=lvl)
    ec_ctx = SQI.ec_ctx

    for k in req:
        if k_test != None and k != k_test:
            continue
        print(f'Level {lvl} - KAT {k}...', end='', flush=True)
        req_k = req[k]
        rsp_k = rsp[k]

        for val in req_k:
            assert req_k[val] == rsp_k[val], f'[Test {k}] {val} not matching'

        seed = req_k["seed"]
        msg = bytes.fromhex(req_k["msg"])
        assert len(msg) == req_k['mlen'], f'[Test {k}] message length {len(msg)} does not match mlen'

        # Compare
        kat_pk = PublicKey.from_hex(rsp_k["pk"], ec_ctx)
        kat_sk = SecretKey.from_hex(rsp_k["sk"], ec_ctx)
        kat_sm = Signature.from_hex(rsp_k["sm"], rsp_k["mlen"], ec_ctx)

        # Sanity on secret key
        assert kat_pk == kat_sk.pk

        # Sanity on signature
        assert kat_sm.msg == msg
        assert sqi.verify(kat_pk, kat_sm)

        # Per-count DRBG seeding, the analogue of the KAT program's
        # randombytes_init(seed, NULL, 256) (apps/PQCgenKAT_sign.c)
        drbg = AES256_CTR_DRBG(bytes.fromhex(seed))
        SQI.randombytes = drbg.random_bytes

        # Keygen (compute pk and sk)
        SQI.keygen()

        # KAT compliance
        assert kat_pk == SQI.pk, "pk mismatch"
        assert kat_sk == SQI.sk, "sk mismatch"

        # Signing
        signature = SQI.sign(msg)
        assert kat_sm == signature

        print(' Done')

    return True


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate SQIsign KATs. Without arguments, all KATs of all levels are checked."
    )
    parser.add_argument("-l", "--level", type=int, choices=[1, 3, 5],
                        help="security level (default: 1)")
    parser.add_argument("-k", "--kat", type=int,
                        help="KAT number (default: all KATs of the given level)")
    args = parser.parse_args()

    if args.level is None and args.kat is None:
        levels = [1, 3, 5]
    else:
        levels = [args.level if args.level is not None else 1]

    for lvl in levels:
        validate(lvl, args.kat)
