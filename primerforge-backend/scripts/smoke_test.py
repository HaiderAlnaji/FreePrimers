#!/usr/bin/env python3
"""
Smoke test for a running PrimerForge backend.

Hits every endpoint with known inputs and checks the response shape
and a few sanity bounds (not full correctness proofs, but enough to
catch "the server is broken" vs "the server works"). Exits non-zero
on any failure so it's CI-friendly.

Usage:
    python3 smoke_test.py http://localhost:8000
    python3 smoke_test.py http://localhost:8000 --skip-specificity
    python3 smoke_test.py http://localhost:8000 --database mirbase_mature --backend local
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request

PASS, FAIL = "\033[32mPASS\033[0m", "\033[31mFAIL\033[0m"


def post(base: str, path: str, body: dict) -> tuple[int, dict | None]:
    req = urllib.request.Request(
        f"{base}{path}",
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            return resp.status, json.loads(resp.read())
    except urllib.error.HTTPError as e:
        try:
            return e.code, json.loads(e.read())
        except Exception:
            return e.code, None
    except Exception as e:
        return 0, {"_exception": str(e)}


def get(base: str, path: str) -> tuple[int, dict | None]:
    try:
        with urllib.request.urlopen(f"{base}{path}", timeout=30) as resp:
            return resp.status, json.loads(resp.read())
    except urllib.error.HTTPError as e:
        try:
            return e.code, json.loads(e.read())
        except Exception:
            return e.code, None
    except Exception as e:
        return 0, {"_exception": str(e)}


def check(label: str, ok: bool, detail: str = "") -> bool:
    print(f"  [{PASS if ok else FAIL}] {label}" + (f"  \u2014 {detail}" if detail and not ok else ""))
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("base_url")
    ap.add_argument("--database", default="mirbase_mature")
    ap.add_argument("--backend", default="local", choices=["local", "ncbi", "auto"])
    ap.add_argument("--skip-specificity", action="store_true")
    args = ap.parse_args()
    base = args.base_url.rstrip("/")

    all_ok = True

    print(f"\nTarget: {base}\n")

    # --- health ---
    print("health")
    status, data = get(base, "/health")
    all_ok &= check("GET /health returns 200", status == 200, f"status={status}")
    all_ok &= check("status field is 'ok'", isinstance(data, dict) and data.get("status") == "ok", str(data))

    # --- thermo/tm ---
    print("\nthermo/tm")
    status, data = post(base, "/thermo/tm", {"seq": "AGCCTAGGATCCGATTATC"})
    ok = status == 200 and isinstance(data, dict) and "tm_c" in data
    all_ok &= check("POST /thermo/tm returns tm_c", ok, json.dumps(data))
    if ok:
        in_range = 30 < data["tm_c"] < 80
        all_ok &= check(f"Tm in plausible range (got {data['tm_c']}\u00b0C)", in_range)

    # --- thermo/hairpin ---
    print("\nthermo/hairpin")
    seq_hairpin = "AGCCTAGGATCCGATTATCGGATCCTAGGCT"  # designed self-complementary
    status, data = post(base, "/thermo/hairpin", {"seq": seq_hairpin})
    ok = status == 200 and isinstance(data, dict) and "dg_kcal_mol" in data
    all_ok &= check("POST /thermo/hairpin returns dg_kcal_mol", ok, json.dumps(data))
    if ok:
        all_ok &= check(f"Structure found on a self-complementary sequence (dG={data['dg_kcal_mol']})",
                         data.get("structure_found") is True and data["dg_kcal_mol"] < 0)

    # --- thermo/duplex ---
    print("\nthermo/duplex")
    status, data = post(base, "/thermo/duplex", {
        "seq1": "AGCCTAGGATCCGATTATC", "seq2": "GATAATCGGATCCTAGGCT",
    })
    ok = status == 200 and isinstance(data, dict) and "dg_kcal_mol" in data
    all_ok &= check("POST /thermo/duplex returns dg_kcal_mol", ok, json.dumps(data))
    if ok:
        all_ok &= check(f"Strong duplex detected on perfect complements (dG={data['dg_kcal_mol']})",
                         data["dg_kcal_mol"] < -10)

    # --- thermo/rna-fold ---
    print("\nthermo/rna-fold")
    status, data = post(base, "/thermo/rna-fold", {"seq": "UAGCUUAUCAGACUGAUGUUGA"})
    ok = status == 200 and isinstance(data, dict) and "structure" in data
    all_ok &= check("POST /thermo/rna-fold returns structure", ok, json.dumps(data))

    # --- thermo/panel-screen ---
    print("\nthermo/panel-screen")
    status, data = post(base, "/thermo/panel-screen", {"seqs": {
        "FI": "AGCCTAGGATCCGATTATC", "RI": "GATAATCGGATCCTAGGCT",
        "FO": "CTGAGTCAGGTTACAGCCTA", "RO": "ACGTTAGCCGTAATCGGA",
    }})
    ok = status == 200 and isinstance(data, list) and len(data) == 10  # C(4,2)+4 = 10
    all_ok &= check(f"POST /thermo/panel-screen returns all 10 pairs (got {len(data) if isinstance(data, list) else 'n/a'})", ok)

    # --- validation ---
    print("\nvalidation")
    status, data = post(base, "/thermo/tm", {"seq": "AGCXYZ"})
    all_ok &= check("Invalid bases rejected with 422", status == 422, f"status={status}")

    # --- specificity ---
    if not args.skip_specificity:
        print(f"\nspecificity (database='{args.database}', backend='{args.backend}')")
        print("  (this can take a while on backend='ncbi' \u2014 NCBI queues can run minutes)")
        status, data = post(base, "/specificity", {
            "seq": "AGCCTAGGATCCGATTATC",
            "database": args.database,
            "backend": args.backend,
            "max_hits": 5,
        })
        if status == 400:
            print(f"  [INFO] Backend reports: {data.get('detail') if data else status}")
            print("  [INFO] This is expected if you haven't run scripts/build_databases.sh yet.")
        else:
            ok = status == 200 and isinstance(data, dict) and "hits" in data
            all_ok &= check("POST /specificity returns a hit list", ok, json.dumps(data)[:300])
            if ok:
                all_ok &= check(f"backend_used reported ({data.get('backend_used')})", "backend_used" in data)
    else:
        print("\nspecificity: skipped (--skip-specificity)")

    print()
    if all_ok:
        print(f"{PASS} \u2014 all checks passed.\n")
        return 0
    else:
        print(f"{FAIL} \u2014 see failures above.\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
