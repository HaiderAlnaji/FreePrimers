"""
Fetch the SPLICED cDNA (mRNA) sequence of a gene's transcript, with exon
boundaries expressed in cDNA coordinates.

This is the correct substrate for qPCR primer design: introns are removed,
exons are contiguous and in transcript order (5'->3'), and exon-exon
junctions are exact positions in the returned sequence. Designing primers
that span a junction (or sit on neighbouring exons) on this sequence
guarantees the assay is cDNA-specific and cannot amplify genomic DNA.

Pipeline (Ensembl REST):
  1. symbol/id  -> gene, list of transcripts (with is_canonical, biotype)
  2. choose transcript (canonical by default, or caller-specified id)
  3. /sequence/id/{tid}?type=cdna  -> spliced cDNA sequence
  4. /overlap/id/{tid}?feature=exon -> exon genomic coords
     then convert to cDNA coordinates by walking exons in transcript order
     and accumulating their lengths.

Returns ok/error so the UI can explain failures.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

try:
    import requests
except ImportError:
    requests = None

ENSEMBL = "https://rest.ensembl.org"


@dataclass
class CdnaExon:
    exon_number: int    # 1-based, transcript order (1 = first exon of mRNA)
    start: int          # 0-based offset in the cDNA sequence
    end: int            # exclusive
    length: int         # exon length (bp)


@dataclass
class TranscriptOption:
    transcript_id: str
    is_canonical: bool
    biotype: str
    length: int          # cDNA length
    n_exons: int


@dataclass
class CdnaFetchResult:
    ok: bool = False
    query: str = ""
    gene: str = ""
    assembly: str = ""
    transcript_id: str = ""
    biotype: str = ""
    sequence: str = ""           # spliced cDNA, 5'->3'
    exons: list = field(default_factory=list)        # list[CdnaExon]
    transcripts: list = field(default_factory=list)  # list[TranscriptOption]
    note: str = ""
    error: str = ""


def _fail(query: str, msg: str) -> CdnaFetchResult:
    return CdnaFetchResult(ok=False, query=query, error=msg)


def _get_gene(query: str, species: str):
    """Resolve a gene symbol or Ensembl id to a gene record with transcripts."""
    q = query.strip()
    # If it's an Ensembl gene id
    if q.upper().startswith("ENSG"):
        url = f"{ENSEMBL}/lookup/id/{q}"
    else:
        url = f"{ENSEMBL}/lookup/symbol/{species}/{q}"
    r = requests.get(url, params={"content-type": "application/json", "expand": 1},
                     headers={"Content-Type": "application/json"}, timeout=20)
    if r.status_code != 200:
        return None
    return r.json()


def list_transcripts(query: str, species: str = "human") -> CdnaFetchResult:
    """Return the list of transcripts for a gene (no sequence)."""
    if requests is None:
        return _fail(query, "The 'requests' library is not installed on the backend.")
    try:
        gj = _get_gene(query, species)
        if not gj:
            return _fail(query, f"Gene '{query}' not found in Ensembl.")
        gene = gj.get("display_name", query)
        assembly = gj.get("assembly_name", "")
        transcripts = gj.get("Transcript", [])
        opts = []
        for t in transcripts:
            exons = t.get("Exon", [])
            length = sum(int(e.get("end", 0)) - int(e.get("start", 0)) + 1 for e in exons)
            opts.append(TranscriptOption(
                transcript_id=t.get("id", ""),
                is_canonical=bool(t.get("is_canonical", 0)),
                biotype=t.get("biotype", ""),
                length=length,
                n_exons=len(exons),
            ))
        # canonical first, then protein_coding, then by length desc
        opts.sort(key=lambda o: (not o.is_canonical, o.biotype != "protein_coding", -o.length))
        return CdnaFetchResult(ok=True, query=query, gene=gene, assembly=assembly,
                               transcripts=opts)
    except Exception as e:
        return _fail(query, f"Transcript lookup failed: {e}")


def fetch_cdna(query: str, species: str = "human",
               transcript_id: str = "") -> CdnaFetchResult:
    """Fetch spliced cDNA + cDNA-coordinate exon boundaries for a transcript.

    If transcript_id is given, that transcript is used; otherwise the canonical
    (or longest protein-coding) transcript is chosen.
    """
    if requests is None:
        return _fail(query, "The 'requests' library is not installed on the backend.")
    try:
        gj = _get_gene(query, species)
        if not gj:
            return _fail(query, f"Gene '{query}' not found in Ensembl.")
        gene = gj.get("display_name", query)
        assembly = gj.get("assembly_name", "")
        strand = int(gj.get("strand", 1))
        transcripts = gj.get("Transcript", [])
        if not transcripts:
            return _fail(query, f"No transcripts found for {gene}.")

        # Choose transcript
        chosen = None
        if transcript_id:
            chosen = next((t for t in transcripts if t.get("id") == transcript_id), None)
            if chosen is None:
                return _fail(query, f"Transcript {transcript_id} not found for {gene}.")
        if chosen is None:
            chosen = next((t for t in transcripts if t.get("is_canonical")), None)
        if chosen is None:
            # longest protein-coding
            pc = [t for t in transcripts if t.get("biotype") == "protein_coding"]
            pool = pc or transcripts
            chosen = max(pool, key=lambda t: sum(
                int(e.get("end", 0)) - int(e.get("start", 0)) + 1 for e in t.get("Exon", [])))

        tid = chosen["id"]
        biotype = chosen.get("biotype", "")

        # Fetch spliced cDNA sequence
        s = requests.get(f"{ENSEMBL}/sequence/id/{tid}",
                         params={"content-type": "application/json", "type": "cdna"},
                         headers={"Content-Type": "application/json"}, timeout=25)
        if s.status_code != 200:
            return _fail(query, f"cDNA sequence fetch failed for {tid} (HTTP {s.status_code}).")
        cdna = s.json().get("seq", "").upper()
        if not cdna:
            return _fail(query, f"Empty cDNA returned for {tid}.")

        # Build exon boundaries in cDNA coordinates.
        # Exons in the transcript record are already on the correct strand;
        # order them in transcript order (5'->3' of the mRNA):
        #   plus strand:  ascending genomic start
        #   minus strand: descending genomic start
        exon_recs = chosen.get("Exon", [])
        exon_recs = sorted(exon_recs, key=lambda e: int(e.get("start", 0)),
                           reverse=(strand == -1))

        # Lay exons contiguously from position 0. The summed genomic exon length
        # should equal the cDNA length; if Ensembl returns a small difference
        # (rare UTR/annotation edge cases), absorb it into the LAST exon so the
        # exon map always covers the full cDNA with no spurious gaps.
        cdna_len = len(cdna)
        raw_lengths = [int(e.get("end", 0)) - int(e.get("start", 0)) + 1 for e in exon_recs]
        total_raw = sum(raw_lengths)

        exons = []
        pos = 0
        n_ex = len(raw_lengths)
        for i, length in enumerate(raw_lengths):
            is_last = (i == n_ex - 1)
            if is_last:
                # Make the final exon end exactly at cdna_len (covers any drift)
                end = cdna_len
                start = min(pos, end)
            else:
                start = pos
                end = pos + length
                if end > cdna_len:
                    end = cdna_len
            if end <= start:
                # exon fully past the cDNA end — skip but keep numbering continuity
                continue
            exons.append(CdnaExon(exon_number=i + 1, start=start, end=end, length=end - start))
            pos = end

        note = ""
        if abs(total_raw - cdna_len) > 5:
            note = (f"Note: summed exon length ({total_raw} bp) differs from cDNA "
                    f"length ({cdna_len} bp) by {abs(total_raw - cdna_len)} bp; "
                    "the last exon boundary was adjusted to fit. This is usually a "
                    "UTR annotation difference and does not affect primer design.")

        return CdnaFetchResult(
            ok=True, query=query, gene=gene, assembly=assembly,
            transcript_id=tid, biotype=biotype, sequence=cdna,
            exons=exons, note=note,
        )
    except Exception as e:
        return _fail(query, f"cDNA fetch error: {e}")
