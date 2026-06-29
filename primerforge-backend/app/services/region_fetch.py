"""
Fetch a genomic region sequence by gene name or by coordinate range (Ensembl).
For the methylation tab: a user gives a gene symbol (e.g. "BRCA2"), an Ensembl
gene id, or a coordinate range ("13:32315000-32316000"), and we return the
reference sequence so they can design bisulfite assays over a promoter/exon.

Network required; returns an object with ok/error so the UI can show why a
lookup failed instead of a generic message.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

try:
    import requests
except ImportError:
    requests = None

ENSEMBL = "https://rest.ensembl.org"


@dataclass
class ExonInfo:
    exon_number: int       # 1-based exon number in transcript
    start: int             # genomic start (0-based offset from region start)
    end: int               # genomic end (exclusive)
    length: int            # exon length in bp

@dataclass
class RegionFetchResult:
    ok: bool = False
    query: str = ""
    sequence: str = ""
    region: str = ""          # resolved "chr:start-end"
    gene: str = ""
    assembly: str = ""
    note: str = ""
    error: str = ""
    exons: list = None        # list of ExonInfo (relative to region start) if available
    transcript_id: str = ""   # Ensembl transcript used for exon lookup
    region_start: int = 0     # absolute genomic start of the returned region


def _fail(q, msg):
    return RegionFetchResult(ok=False, query=q, error=msg)


def _looks_like_region(q: str) -> bool:
    # e.g. 13:32315000-32316000  or  13:32315000..32316000
    import re
    return bool(re.match(r"^[\dXYMT]+:\d+[-.]+\d+$", q.strip()))


def _fetch_exons(species: str, gene_id: str, region_start: int, region_end: int) -> tuple:
    """Fetch exon coordinates for the canonical transcript of a gene.
    Returns (transcript_id, [ExonInfo]) with positions relative to region_start.
    ExonInfo.exon_number is the transcript exon number (1 = first exon of mRNA).
    Returns ("", []) on any failure (non-fatal).
    """
    try:
        # Get canonical transcript and gene strand
        r = requests.get(
            f"{ENSEMBL}/lookup/id/{gene_id}",
            params={"content-type": "application/json", "expand": 1},
            headers={"Content-Type": "application/json"}, timeout=20)
        if r.status_code != 200:
            return "", []
        gj = r.json()
        gene_strand = int(gj.get("strand", 1))   # +1 or -1
        transcripts = gj.get("Transcript", [])
        canon = next((t for t in transcripts if t.get("is_canonical")), None)
        if not canon and transcripts:
            canon = max(transcripts, key=lambda t: t.get("Translation", {}).get("length", 0) if t.get("Translation") else 0)
        if not canon:
            return "", []
        tid = canon["id"]

        # Get exons for this transcript
        er = requests.get(
            f"{ENSEMBL}/overlap/id/{tid}",
            params={"content-type": "application/json", "feature": "exon"},
            headers={"Content-Type": "application/json"}, timeout=20)
        if er.status_code != 200:
            return tid, []
        exon_data = er.json()

        # Sort in TRANSCRIPT order:
        # Plus strand  (+1): ascending genomic start  → exon 1 is leftmost
        # Minus strand (-1): descending genomic start → exon 1 is rightmost
        exon_data = sorted(exon_data, key=lambda e: e.get("start", 0),
                           reverse=(gene_strand == -1))

        exons = []
        for i, e in enumerate(exon_data):
            es = int(e.get("start", 0))
            ee = int(e.get("end", 0))
            # Convert to 0-based offset from region_start
            # For minus-strand genes Ensembl returns the reverse-complement of the
            # genomic region, so position in the returned sequence for a minus-strand
            # gene is: seq_pos = (region_end - genomic_end) for the exon start.
            if gene_strand == -1:
                # The fetched sequence is the reverse complement of [region_start, region_end]
                # Exon genomic [es, ee] maps to sequence positions:
                #   seq_start = region_end - ee   (0-based from left of returned seq)
                #   seq_end   = region_end - es + 1
                seq_start = region_end - ee
                seq_end   = region_end - es + 1
            else:
                seq_start = es - region_start
                seq_end   = ee - region_start + 1

            # Clip to fetched region
            rel_start = max(0, seq_start)
            rel_end   = min(region_end - region_start, seq_end)
            if rel_start >= rel_end:
                continue
            exons.append(ExonInfo(
                exon_number=i + 1,   # 1-based transcript exon number
                start=rel_start,
                end=rel_end,
                length=ee - es + 1,
            ))
        return tid, exons
    except Exception:
        return "", []


def fetch_region(query: str, species: str = "human", flank: int = 0,
                 max_len: int = 4000, fetch_exons: bool = True) -> RegionFetchResult:
    query = (query or "").strip()
    if not query:
        return _fail(query, "No gene or region given.")
    if requests is None:
        return _fail(query, "The 'requests' library is not installed on the backend.")

    region = None
    gene_name = ""
    assembly = ""
    gene_id_for_exons = ""

    if _looks_like_region(query):
        region = query.replace("..", "-")
    else:
        # treat as a gene symbol or Ensembl id -> look up its coordinates
        try:
            g = requests.get(f"{ENSEMBL}/lookup/symbol/{species}/{query}",
                             params={"content-type": "application/json"},
                             headers={"Content-Type": "application/json"}, timeout=20)
            if g.status_code == 404:
                # maybe it's an Ensembl stable id
                g = requests.get(f"{ENSEMBL}/lookup/id/{query}",
                                 params={"content-type": "application/json"},
                                 headers={"Content-Type": "application/json"}, timeout=20)
            if g.status_code != 200:
                return _fail(query, f"Gene '{query}' not found in Ensembl (HTTP {g.status_code}).")
            gj = g.json()
            chrom = gj.get("seq_region_name")
            start = int(gj.get("start"))
            end = int(gj.get("end"))
            gene_name = gj.get("display_name", query)
            assembly = gj.get("assembly_name", "")
            gene_id_for_exons = gj.get("id", "")
            strand = int(gj.get("strand", 1))
            span = min(max_len, end - start + 1)
            if strand == 1:
                rs, re_ = start - flank, start + span - 1
            else:
                rs, re_ = end - span + 1, end + flank
            region = f"{chrom}:{rs}-{re_}"
        except Exception as e:
            return _fail(query, f"Gene lookup failed: {e}")

    # parse region start for exon offset calculation
    try:
        _, coords = region.split(":")
        r_start, r_end = [int(x) for x in coords.replace("..", "-").split("-")]
    except Exception:
        r_start, r_end = 0, 0

    # fetch the sequence for the resolved region
    try:
        s = requests.get(f"{ENSEMBL}/sequence/region/{species}/{region}",
                         params={"content-type": "application/json"},
                         headers={"Content-Type": "application/json"}, timeout=25)
        if s.status_code != 200:
            return _fail(query, f"Sequence fetch failed for {region} (HTTP {s.status_code}).")
        seq = s.json().get("seq", "").upper()
        if not seq:
            return _fail(query, f"Empty sequence returned for {region}.")
        if len(seq) > max_len:
            seq = seq[:max_len]
            note = f"Sequence truncated to {max_len} bp."
        else:
            note = ""
        # fetch exon positions if we have a gene id
        exons, tid = [], ""
        if fetch_exons and gene_id_for_exons:
            tid, exons = _fetch_exons(species, gene_id_for_exons, r_start, r_end)
        return RegionFetchResult(ok=True, query=query, sequence=seq, region=region,
                                 gene=gene_name, assembly=assembly, note=note,
                                 exons=exons if exons else None,
                                 transcript_id=tid, region_start=r_start)
    except Exception as e:
        return _fail(query, f"Sequence fetch error: {e}")
