"""
SNP-aware primer placement scorer.

A primer that overlaps a common population variant can fail to bind in carriers
of the alternate allele, causing allele dropout (failed amplification / biased
genotype) or a Cq shift in qPCR. This module scores any candidate primer for
the variant-overlap risk so a designer (standard PCR/qPCR OR the ARMS/miRNA
engines) can avoid placing primers -- especially their 3' ends -- over common
SNPs.

Empirical basis for the position weighting (not a guessed curve):
  - Mismatches within the last ~5 bases of the 3' end have the most dramatic
    effect; the 3'-terminal base shifts qPCR Cq by up to ~7 (128-fold)
    (Lefever et al., via IDT 2024 technical note).
  - Mismatches in the last 3-4 positions give minimal/no extension; INTERNAL
    mismatches cause preferential amplification (allele bias) rather than total
    dropout -- e.g. a documented BRCA1 false-negative from a SNP at the 11th
    base from the 3' end (Silva et al. 2017, Mol Genet Genomic Med).
  - Purine/purine and pyrimidine/pyrimidine 3'-terminal mismatches (esp. A/G,
    C/C) are the most disruptive base compositions (Lefever et al.).

Scoring (per overlapping variant):
    risk = base_position_weight(distance_from_3prime) * frequency_factor(MAF)
The primer's total risk is the sum over overlapping variants; a single
high-frequency 3'-proximal variant dominates, as it should.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


@dataclass
class Variant:
    pos: int            # 0-based position in the SAME coordinate frame as the region
    ref: str
    alt: str
    maf: float = 0.0    # minor/alternate allele frequency (0..0.5+); 0 if unknown
    rsid: str = ""


@dataclass
class PrimerRisk:
    seq: str
    start: int          # 0-based start in region coords
    end: int            # exclusive
    strand: str         # 'fwd' or 'rev'
    overlaps: list = field(default_factory=list)   # (Variant, dist_from_3prime, risk)
    total_risk: float = 0.0
    verdict: str = ""   # 'clean' | 'caution' | 'avoid'


# --- position weighting (distance measured from the primer 3' end, 0 = terminus) ---
def position_weight(dist_from_3prime: int) -> float:
    """Relative disruption of a mismatch at a given distance from the 3' end.

    Anchored to the empirical picture: terminal base ~1.0 (max), steep through
    the last 5 bases, a smaller plateau internally (bias, not dropout), and
    near-zero far into the 5' region.
    """
    if dist_from_3prime < 0:
        return 0.0
    if dist_from_3prime == 0:
        return 1.0          # 3'-terminal: strongest (up to ~7 Cq)
    if dist_from_3prime == 1:
        return 0.85         # penultimate
    if dist_from_3prime == 2:
        return 0.65
    if dist_from_3prime <= 4:
        return 0.45         # within last 5 bases: still dramatic
    if dist_from_3prime <= 9:
        return 0.20         # internal: preferential amplification / bias
    return 0.08             # deep 5' side: usually tolerated


def frequency_factor(maf: float) -> float:
    """Scale risk by how often the alternate allele is actually present.

    A variant nobody carries is harmless; a common one is dangerous. Unknown
    MAF (0.0) is treated as a moderate default so unannotated variants are not
    silently ignored.
    """
def frequency_factor(maf: float, unknown_weight: float = 0.15) -> float:
    """Scale risk by how often the alternate allele is actually present.

    A variant nobody carries is harmless; a common one is dangerous.

    unknown_weight: weight applied when MAF is 0.0/unknown. Two regimes:
      - PLACEMENT (default, 0.15): for routine PCR/qPCR in a general
        population, only KNOWN-COMMON variants threaten amplification; a
        cataloged-but-rare/unannotated variant (e.g. a clinical BRCA2 indel
        seen once) essentially never causes dropout in a random sample, so it
        is treated as low risk. This avoids the failure mode where a heavily
        clinically-curated gene looks "unusable" purely because every rare
        mutation ever seen is logged there.
      - CONSERVATIVE (0.5, set by caller): for clinical/diagnostic design where
        you want to flag ANY cataloged variant under the primer regardless of
        frequency.
    Note: a TRUE common variant with MAF=0 would be mis-handled, but by
    definition common variants have measured MAF; MAF=0 means rare or
    unannotated, which is exactly what we down-weight for placement.
    """
    if maf <= 0.0:
        return unknown_weight   # unknown/unannotated -> caller-chosen (default low)
    if maf < 0.001:
        return 0.15             # very rare
    if maf < 0.01:
        return 0.4
    if maf < 0.05:
        return 0.7
    return 1.0                  # common (>=5%): full weight


# terminal base-composition multiplier (most disruptive 3'-terminal mispairs)
def _terminal_composition_multiplier(primer_3p_base: str, alt_base_on_template: str) -> float:
    """If the variant sits exactly at the 3' terminus, weight by mismatch type.

    A/G (purine/purine) and C/C (pyrimidine/pyrimidine) are the most disruptive
    (largest Cq shift); transitions are milder. Returns ~0.8..1.2.
    """
    pair = {primer_3p_base, alt_base_on_template}
    if pair == {"A", "G"} or pair == {"C"}:
        return 1.2
    if pair == {"C", "T"}:
        return 0.85          # C.T extends relatively readily
    return 1.0


def score_primer(seq: str, start: int, strand: str, variants: list,
                 region_offset: int = 0, unknown_weight: float = 0.15) -> PrimerRisk:
    """Score one primer for variant-overlap risk.

    seq      : primer sequence as it would be ordered (5'->3').
    start    : 0-based start of the primer's footprint on the region (the
               low-coordinate end on the sense strand), in region coords.
    strand   : 'fwd' (3' end at high coord) or 'rev' (3' end at low coord).
    variants : list[Variant] in region coords (after region_offset applied).
    """
    L = len(seq)
    end = start + L
    pr = PrimerRisk(seq=seq, start=start, end=end, strand=strand)
    for v in variants:
        vp = v.pos - region_offset
        if not (start <= vp < end):
            continue
        # distance of this variant from the primer's 3' end
        if strand == "fwd":
            dist_3p = (end - 1) - vp          # 3' end is the high-coord base
            primer_3p_base = seq[-1]
        else:
            dist_3p = vp - start              # 3' end is the low-coord base (rev)
            primer_3p_base = seq[-1]
        w = position_weight(dist_3p)
        f = frequency_factor(v.maf, unknown_weight=unknown_weight)
        risk = w * f
        if dist_3p == 0:
            risk *= _terminal_composition_multiplier(primer_3p_base, v.alt)
        pr.overlaps.append((v, dist_3p, round(risk, 3)))

    # --- aggregation: a primer fails mainly because of its WORST single variant,
    # not because many independent rare variants happen to be cataloged under it.
    # Naive summation scales with database curation density (a heavily studied
    # gene like BRCA2 would look "unusable" purely from logged rare mutations),
    # which is wrong. Use the dominant variant plus a damped contribution from
    # the rest: total = max + 0.25 * sum(others). This keeps a single common /
    # 3'-proximal variant decisive while not letting a pile of MAF~0 variants
    # accumulate into a false "avoid".
    contribs = sorted((r for (_, _, r) in pr.overlaps), reverse=True)
    if contribs:
        pr.total_risk = round(contribs[0] + 0.25 * sum(contribs[1:]), 3)
    else:
        pr.total_risk = 0.0
    # verdict driven by the WORST single variant's risk and its position
    worst_single = contribs[0] if contribs else 0.0
    worst_dist = min((d for (_, d, _) in pr.overlaps), default=99)
    if not pr.overlaps:
        pr.verdict = "clean"
    elif worst_single >= 0.5 or (worst_dist <= 1 and worst_single >= 0.3):
        pr.verdict = "avoid"          # a genuinely dangerous single variant
    elif worst_single >= 0.2 or pr.total_risk >= 0.6:
        pr.verdict = "caution"
    else:
        pr.verdict = "clean"
    return pr


def scan_region(region_seq: str, variants: list, primer_len: int = 20,
                strand: str = "fwd", region_offset: int = 0,
                top: int = 10, unknown_weight: float = 0.15) -> list:
    """Slide a primer-sized window across the region and rank windows by lowest
    variant-overlap risk (best primer-placement positions).
    Returns the `top` lowest-risk PrimerRisk objects.
    """
    region_seq = region_seq.upper()
    out = []
    for s in range(0, len(region_seq) - primer_len + 1):
        sub = region_seq[s:s + primer_len]
        if strand == "rev":
            ordered = "".join(COMP[b] for b in reversed(sub))
        else:
            ordered = sub
        pr = score_primer(ordered, s, strand, variants, region_offset, unknown_weight)
        out.append(pr)
    out.sort(key=lambda p: (p.total_risk, p.start))
    return out[:top]


def risk_report(pr: PrimerRisk) -> str:
    tag = {"clean": "CLEAN", "caution": "CAUTION", "avoid": "AVOID"}.get(pr.verdict, "?")
    lines = [f"[{tag}] primer 5'-{pr.seq}-3' ({pr.strand}, region {pr.start}-{pr.end}) "
             f"total risk {pr.total_risk}"]
    for (v, d, r) in sorted(pr.overlaps, key=lambda x: x[1]):
        where = "3'-TERMINUS" if d == 0 else f"{d} nt from 3'"
        rs = f" {v.rsid}" if v.rsid else ""
        lines.append(f"    variant{rs} {v.ref}>{v.alt} MAF={v.maf} at {where} -> risk {r}")
    return "\n".join(lines)
