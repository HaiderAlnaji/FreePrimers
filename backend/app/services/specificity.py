"""
Specificity service.

Checks where a primer/probe sequence binds across a reference set
(miRBase mature miRNAs, a transcriptome, a genome region, etc.) using
BLAST. Two backends:

- local: blastn against a pre-built local BLAST database. Fast,
  unlimited, works offline, and is what you want for any production
  or batch use. Requires `makeblastdb`-built databases under
  settings.blast_db_dir (see scripts/build_databases.sh).
- ncbi: the NCBI Common URL API (Put/Get). No local setup, but slow
  (each search queues on NCBI's servers, typically 30s-several
  minutes), rate-limited, and intended for occasional/interactive use
  only — do not batch against it.

backend="auto" tries local first and falls back to NCBI only if no
local database with that name is configured.
"""

from __future__ import annotations

import asyncio
import csv
import io
import re
import subprocess
import time
from pathlib import Path
from typing import Optional

import httpx

from app.config import settings
from app.models.schemas import SpecificityHit, SpecificityRequest, SpecificityResult

# Short oligos (<30 nt) need blastn-short / adjusted word size, or
# nothing will hit. This mirrors NCBI's own short_query_adjust logic.
SHORT_QUERY_THRESHOLD = 30

# Tabular output columns requested from blastn -outfmt, and from the
# NCBI URL API's FORMAT_TYPE=Tabular / ALIGNMENT_VIEW=Tabular output.
# Keeping these aligned means both backends parse through one function.
TABULAR_FIELDS = [
    "sseqid", "stitle", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


class SpecificityError(Exception):
    pass


# Accession prefixes that are definitively miRNA/ncRNA records
_MIRNA_ACC_PREFIXES = ("MI", "MIMAT", "XR_", "NR_")

# Keywords in subject titles that identify miRNA/ncRNA records
_MIRNA_TITLE_KEYWORDS = (
    "microrna", "mirna", "mir-", "let-7", "hsa-", "mmu-", "rno-",
    "small rna", "ncrna", "non-coding", "noncoding",
)

# Mature miRNAs are 18–25 nt; forward primers are the miRNA + optional 5' tail (≤28 nt).
# The stem-loop RT primer (≥45 nt) and universal reverse are longer — apply standard
# risk logic to those so genuine off-targets in those primers are still flagged.
_MIRNA_QUERY_LEN_THRESHOLD = 28


def _is_mirna_record(subject_id: str, subject_title: str) -> bool:
    """Return True if the hit looks definitively like a miRNA/ncRNA record."""
    sid = (subject_id or "").upper()
    stitle = (subject_title or "").lower()
    for pfx in _MIRNA_ACC_PREFIXES:
        if sid.startswith(pfx.upper()):
            return True
    for kw in _MIRNA_TITLE_KEYWORDS:
        if kw in stitle:
            return True
    return False


def _classify_risk(pident: float, align_length: int, query_len: int, qend: int,
                   is_mirna: bool = False,
                   subject_id: str = "", subject_title: str = "") -> tuple[str, str]:
    """Return (risk_level, hit_context).

    Standard primer logic: high/moderate/low based on identity, coverage, 3′ reach.

    miRNA primer logic (is_mirna=True):
    - If the accession/title identifies a miRNA/ncRNA record → expected.
    - If the query is short (≤_MIRNA_QUERY_LEN_THRESHOLD) AND the hit covers ≥90%
      of the query → expected.  A mature miRNA or stem-loop RT primer (18–50 nt)
      that matches at full coverage is hitting the miRNA itself or a family member
      regardless of the accession type; this is not off-target amplification.
    - Otherwise: apply standard risk logic.  A partial hit to a long mRNA or a
      full hit to an unrelated coding sequence is a real concern.
    """
    coverage = align_length / query_len if query_len else 0
    reaches_3prime = qend >= query_len - 2

    if is_mirna:
        # Definitively-identified miRNA/ncRNA record
        if _is_mirna_record(subject_id, subject_title):
            return "expected", "mirna_family"
        # Short query (mature miRNA / stem-loop primer) with full-length coverage
        # → the hit is the miRNA sequence embedded in whatever transcript record;
        # it is the TARGET, not off-target.
        if query_len <= _MIRNA_QUERY_LEN_THRESHOLD and coverage >= 0.9:
            return "expected", "mirna_family"

    # Standard off-target risk for all other cases
    if coverage >= 0.9 and pident >= 90 and reaches_3prime:
        return "high", "off_target"
    if coverage >= 0.7 and pident >= 80:
        return "moderate", "off_target"
    return "low", "off_target"


def _parse_tabular(raw: str, query_len: int, is_mirna: bool = False) -> list[SpecificityHit]:
    hits: list[SpecificityHit] = []
    reader = csv.reader(io.StringIO(raw.strip()), delimiter="\t")
    for row in reader:
        if not row or row[0].startswith("#"):
            continue
        if len(row) < len(TABULAR_FIELDS):
            continue
        sseqid, stitle, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = row[: len(TABULAR_FIELDS)]
        qend_i = int(qend)
        risk, hit_context = _classify_risk(
            float(pident), int(length), query_len, qend_i,
            is_mirna=is_mirna, subject_id=sseqid, subject_title=stitle or "")
        hit = SpecificityHit(
            subject_id=sseqid,
            subject_title=stitle or None,
            pct_identity=float(pident),
            align_length=int(length),
            mismatches=int(mismatch),
            gaps=int(gapopen),
            query_start=int(qstart),
            query_end=qend_i,
            subject_start=int(sstart),
            subject_end=int(send),
            evalue=float(evalue),
            bit_score=float(bitscore),
            risk=risk,
            hit_context=hit_context,
        )
        hits.append(hit)
    return hits


# ---------------------------------------------------------------------------
# Local BLAST+ backend
# ---------------------------------------------------------------------------

def local_database_path(name: str) -> Path:
    return Path(settings.blast_db_dir) / name


def local_database_exists(name: str) -> bool:
    db = local_database_path(name)
    # nucleotide DBs built with makeblastdb leave a .nhr/.nin/.nsq
    # trio, or a single .ndb for newer BLAST+ versions.
    return db.with_suffix(".nhr").exists() or db.with_suffix(".ndb").exists()


def list_local_databases() -> list[str]:
    d = Path(settings.blast_db_dir)
    if not d.exists():
        return []
    names = set()
    for f in d.iterdir():
        if f.suffix in (".nhr", ".ndb"):
            names.add(f.stem)
    return sorted(names)


async def run_local_blast(req: SpecificityRequest) -> SpecificityResult:
    db_path = local_database_path(req.database)
    if not local_database_exists(req.database):
        raise SpecificityError(
            f"No local BLAST database named '{req.database}' found in "
            f"{settings.blast_db_dir}. Build it with scripts/build_databases.sh "
            f"or use backend='ncbi' / 'auto'."
        )

    word_size = "7" if len(req.seq) < SHORT_QUERY_THRESHOLD else "11"
    cmd = [
        settings.blastn_path,
        "-query", "-",
        "-db", str(db_path),
        "-task", "blastn-short" if len(req.seq) < SHORT_QUERY_THRESHOLD else "blastn",
        "-word_size", word_size,
        "-evalue", "1000" if len(req.seq) < SHORT_QUERY_THRESHOLD else "10",
        "-outfmt", "6 " + " ".join(TABULAR_FIELDS),
        "-max_target_seqs", str(req.max_hits),
    ]
    fasta = f">query\n{req.seq}\n"

    proc = await asyncio.create_subprocess_exec(
        *cmd,
        stdin=asyncio.subprocess.PIPE,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    stdout, stderr = await proc.communicate(fasta.encode())
    if proc.returncode != 0:
        raise SpecificityError(f"blastn failed: {stderr.decode(errors='replace')}")

    hits = _parse_tabular(stdout.decode(), query_len=len(req.seq), is_mirna=getattr(req,"is_mirna",False))
    return SpecificityResult(
        seq=req.seq,
        database=req.database,
        backend_used="local",
        total_hits=len(hits),
        hits=hits[: req.max_hits],
    )


# ---------------------------------------------------------------------------
# NCBI Common URL API backend (fallback)
# ---------------------------------------------------------------------------

NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
NCBI_POLL_INTERVAL_S = 15
NCBI_MAX_WAIT_S = 300  # 5 minutes; short oligo searches are usually faster


async def run_ncbi_blast(req: SpecificityRequest) -> SpecificityResult:
    """
    Submits via CMD=Put, polls CMD=Get until the job leaves WAITING,
    then re-requests with FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular for
    a parseable tab-separated report. This is NCBI's shared, public,
    rate-limited service — keep usage occasional and never call this
    in a tight batch loop (NCBI will rate-limit or block the caller).
    """
    warnings = [
        "Using the public NCBI BLAST API: this can take 30s\u2013several "
        "minutes and is rate-limited. Configure a local database for "
        "production or batch use."
    ]
    async with httpx.AsyncClient(timeout=60.0) as client:
        put_resp = await client.get(
            NCBI_BLAST_URL,
            params={
                "CMD": "Put",
                "PROGRAM": "blastn",
                "MEGABLAST": "off",  # megablast is tuned for long, near-identical sequences
                "DATABASE": req.database,
                "QUERY": req.seq,
                "SHORT_QUERY_ADJUST": "true" if len(req.seq) < SHORT_QUERY_THRESHOLD else "false",
                "HITLIST_SIZE": str(req.max_hits),
            },
        )
        put_resp.raise_for_status()
        rid_match = re.search(r"RID = (\S+)", put_resp.text)
        rtoe_match = re.search(r"RTOE = (\d+)", put_resp.text)
        if not rid_match:
            raise SpecificityError("NCBI did not return a request ID (RID). Response may indicate a malformed query.")
        rid = rid_match.group(1)
        initial_wait = int(rtoe_match.group(1)) if rtoe_match else NCBI_POLL_INTERVAL_S

        await asyncio.sleep(min(initial_wait, NCBI_POLL_INTERVAL_S))
        elapsed = 0
        while elapsed < NCBI_MAX_WAIT_S:
            status_resp = await client.get(
                NCBI_BLAST_URL, params={"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}
            )
            status_resp.raise_for_status()
            if "Status=READY" in status_resp.text:
                if "ThereAreHits=yes" not in status_resp.text:
                    return SpecificityResult(
                        seq=req.seq, database=req.database, backend_used="ncbi",
                        total_hits=0, hits=[], warnings=warnings,
                    )
                break
            if "Status=FAILED" in status_resp.text or "Status=UNKNOWN" in status_resp.text:
                raise SpecificityError(f"NCBI BLAST job {rid} failed or expired.")
            await asyncio.sleep(NCBI_POLL_INTERVAL_S)
            elapsed += NCBI_POLL_INTERVAL_S
        else:
            raise SpecificityError(f"NCBI BLAST job {rid} did not finish within {NCBI_MAX_WAIT_S}s.")

        result_resp = await client.get(
            NCBI_BLAST_URL,
            params={
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "Text",
                "ALIGNMENT_VIEW": "Tabular",
                "DESCRIPTIONS": str(req.max_hits),
                "ALIGNMENTS": str(req.max_hits),
            },
        )
        result_resp.raise_for_status()

    hits = _parse_ncbi_tabular_text(result_resp.text, query_len=len(req.seq),
                                    is_mirna=getattr(req, "is_mirna", False))
    return SpecificityResult(
        seq=req.seq, database=req.database, backend_used="ncbi",
        total_hits=len(hits), hits=hits[: req.max_hits], warnings=warnings,
    )


def _parse_ncbi_tabular_text(text: str, query_len: int, is_mirna: bool = False) -> list[SpecificityHit]:
    """
    NCBI's Text+Tabular report embeds a tab-separated hit table
    inside an HTML/text wrapper. Pull out lines that look like the
    12-column tabular format (same columns as local blastn -outfmt 6,
    NCBI just doesn't add a clean delimiter around the block).
    """
    candidate_lines = [
        line for line in text.splitlines()
        if line.count("\t") >= len(TABULAR_FIELDS) - 1
    ]
    if not candidate_lines:
        return []
    return _parse_tabular("\n".join(candidate_lines), query_len=query_len, is_mirna=is_mirna)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

async def check_specificity(req: SpecificityRequest) -> SpecificityResult:
    if req.backend == "local":
        return await run_local_blast(req)
    if req.backend == "ncbi":
        return await run_ncbi_blast(req)

    # auto: prefer local, fall back to NCBI only on missing-DB, not on
    # e.g. a transient blastn crash (that should surface as an error).
    if local_database_exists(req.database):
        return await run_local_blast(req)
    result = await run_ncbi_blast(req)
    result.warnings = [
        f"No local database '{req.database}' configured; used the NCBI "
        f"API instead.",
        *result.warnings,
    ]
    return result
