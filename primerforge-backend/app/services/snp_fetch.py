"""
Fetch SNP flanking sequence by rsID (Ensembl REST), for the 'reference SNP'
input box. Network access is required and may be unavailable (sandbox, offline,
API down); callers must handle failure and fall back to paste/upload.

Returns a SnpFetchResult always; on failure, ok=False and `error` explains why,
so the UI can show a specific reason instead of a generic "lookup failed".
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

try:
    import requests
except ImportError:
    requests = None

ENSEMBL = "https://rest.ensembl.org"
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


@dataclass
class SnpFetchResult:
    rsid: str
    ok: bool = False
    sense: str = ""
    snp_index: int = 0
    allele1: str = ""
    allele2: str = ""
    assembly: str = ""
    chromosome: str = ""
    position: int = 0
    gene: str = ""
    consequence: str = ""
    maf: float = 0.0
    ancestral: str = ""
    var_class: str = ""
    note: str = ""
    error: str = ""


def _fail(rsid, msg):
    return SnpFetchResult(rsid=rsid, ok=False, error=msg)


def fetch_snp(rsid: str, flank: int = 250, species: str = "human") -> SnpFetchResult:
    rsid = (rsid or "").strip()
    if not rsid:
        return _fail(rsid, "No rsID given.")
    if requests is None:
        return _fail(rsid, "The 'requests' library is not installed on the backend "
                           "(pip install requests).")
    # 1) variation endpoint: location + alleles
    try:
        v = requests.get(
            f"{ENSEMBL}/variation/{species}/{rsid}",
            params={"content-type": "application/json"},
            headers={"Content-Type": "application/json"}, timeout=20)
    except Exception as e:
        return _fail(rsid, f"Could not reach Ensembl (network error): {e}")
    if v.status_code == 404:
        return _fail(rsid, f"'{rsid}' not found in Ensembl (404). Check the rsID.")
    if v.status_code != 200:
        return _fail(rsid, f"Ensembl variation endpoint returned HTTP {v.status_code}.")
    try:
        vj = v.json()
    except Exception:
        return _fail(rsid, "Ensembl returned a non-JSON response for the variation lookup.")
    mappings = vj.get("mappings") or []
    if not mappings:
        return _fail(rsid, f"'{rsid}' has no genomic mappings in Ensembl.")
    m = mappings[0]
    allele_str = m.get("allele_string", "")
    parts = allele_str.split("/")
    if len(parts) < 2:
        return _fail(rsid, f"Unexpected allele string from Ensembl: '{allele_str}'.")
    ref_allele, alt_allele = parts[0], parts[1]
    seq_region = m.get("seq_region_name")
    try:
        start = int(m.get("start"))
    except Exception:
        return _fail(rsid, "Ensembl mapping had no usable start coordinate.")
    strand = int(m.get("strand", 1))
    assembly = m.get("assembly_name", "")
    # extra details for display
    consequence = vj.get("most_severe_consequence", "")
    maf = float(vj.get("MAF") or 0.0)
    ancestral = vj.get("ancestral_allele", "") or ""
    var_class = vj.get("var_class", "") or ""   # e.g. 'SNP', 'insertion', 'deletion'
    # Gene symbol: the /variation endpoint does NOT return overlapping genes, so a
    # second call to the overlap endpoint is needed. Wrapped so a failure here
    # never breaks the (already successful) variant+sequence fetch — gene just
    # stays blank. [VERIFY against live Ensembl: field is 'external_name' on the
    # gene feature; older releases used 'external_name', current returns it too.]
    gene = ""
    try:
        ov = requests.get(
            f"{ENSEMBL}/overlap/region/{species}/{seq_region}:{start}-{start}",
            params={"feature": "gene", "content-type": "application/json"},
            headers={"Content-Type": "application/json"}, timeout=15)
        if ov.status_code == 200:
            genes = ov.json()
            if isinstance(genes, list) and genes:
                # prefer a protein-coding gene if several overlap
                coding = [g for g in genes if g.get("biotype") == "protein_coding"]
                pick = (coding or genes)[0]
                gene = pick.get("external_name") or pick.get("id") or ""
    except Exception:
        gene = ""

    # 2) sequence endpoint: flanking reference
    lo, hi = start - flank, start + flank
    region = f"{seq_region}:{lo}-{hi}:{strand}"
    try:
        s = requests.get(
            f"{ENSEMBL}/sequence/region/{species}/{region}",
            params={"content-type": "application/json"},
            headers={"Content-Type": "application/json"}, timeout=20)
    except Exception as e:
        return _fail(rsid, f"Could not reach Ensembl sequence endpoint: {e}")
    if s.status_code != 200:
        return _fail(rsid, f"Ensembl sequence endpoint returned HTTP {s.status_code}.")
    try:
        seq = s.json().get("seq", "").upper()
    except Exception:
        return _fail(rsid, "Ensembl returned a non-JSON response for the sequence lookup.")
    if not seq:
        return _fail(rsid, "Ensembl returned an empty sequence for the region.")

    snp_index = flank if flank < len(seq) else len(seq) // 2
    a1, a2 = ref_allele, alt_allele
    note = ""
    if snp_index >= len(seq):
        return _fail(rsid, "Returned sequence shorter than expected; cannot locate variant.")
    if seq[snp_index] not in (a1, a2):
        a1c, a2c = COMP.get(a1, a1), COMP.get(a2, a2)
        if seq[snp_index] in (a1c, a2c):
            a1, a2 = a1c, a2c
        else:
            note = (f"base at variant position ({seq[snp_index]}) does not match reported "
                    f"alleles {ref_allele}/{alt_allele}; verify strand before ordering.")
    return SnpFetchResult(rsid=rsid, ok=True, sense=seq, snp_index=snp_index,
                          allele1=a1, allele2=a2, assembly=assembly, note=note,
                          chromosome=seq_region, position=start, gene=gene,
                          consequence=consequence, maf=maf, ancestral=ancestral, var_class=var_class)
