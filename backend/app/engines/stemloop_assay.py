"""
Full stem-loop RT-qPCR assay designer for miRNA family discrimination.

Architecture (Chen et al. 2005; Varkonyi-Gasic & Hellens 2011):
  - STEM-LOOP RT PRIMER = universal backbone (forms hairpin + carries the
    universal reverse-primer site) + a 6-nt 3' extension that is the reverse
    complement of the miRNA's last 6 nt. RT is primed off that 6-nt anneal.
  - FORWARD PRIMER = miRNA-specific, covering the miRNA sequence EXCLUDING the
    last 6 nt (those are occupied by the RT primer), usually with a 5' tail to
    raise Tm.
  - UNIVERSAL REVERSE PRIMER = complementary to a site in the RT backbone;
    identical for every assay.

CRITICAL geometry the discrimination engine must respect:
  The forward primer can only place its 3' terminus on miRNA positions that are
  NOT inside the last-6-nt RT-anneal window. So a discriminating base at:
    - position <= (len - 6): reachable by the forward-primer 3' end -> strong
      kinetic discrimination (problem #1 model) is AVAILABLE.
    - position  > (len - 6): falls in the RT 6-nt window -> discrimination must
      come from the RT step (RNA:DNA), where it is INTERNAL to the 6-mer and
      weak (G.U wobble etc.). The kinetic 3'-block is NOT available unless the
      RT extension is shortened (a real specificity/Tm tradeoff this tool flags
      but does not silently assume).

This module designs all three primers and reports, per sibling, which
discrimination regime applies and the predicted strength.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from . import empirical_discrimination as emp

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Canonical universal stem-loop backbone (Varkonyi-Gasic & Hellens style).
# Lowercase-free; this is the portion 5' of the 6-nt miRNA-specific extension.
# The universal reverse primer anneals within this backbone.
UNIVERSAL_BACKBONE = "GTCGTATCCAGTGCAGGGTCCGAGGTATTCGCACTGGATACGAC"
UNIVERSAL_REVERSE = "GTGCAGGGTCCGAGGT"   # anneals within backbone (a common choice)
RT_EXTENSION_NT = 6


def to_dna(seq: str) -> str:
    return seq.strip().upper().replace("U", "T")


def revcomp(seq: str) -> str:
    return "".join(COMP[b] for b in reversed(to_dna(seq)))


def tm_wallace(seq: str) -> float:
    """Quick Tm estimate (Wallace rule) for short primer Tm sanity only."""
    seq = to_dna(seq)
    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    return 2 * at + 4 * gc


@dataclass
class AssaySibling:
    name: str
    diff_pos: Optional[int]          # 0-based nearest-3' difference
    from_3p: Optional[int]
    regime: str                      # 'forward-terminal' | 'rt-window' | 'internal' | 'identical'
    terminal_mismatch: Optional[str]
    discrimination_log: float
    note: str = ""


@dataclass
class StemLoopAssay:
    target: str
    target_dna: str
    rt_primer: str
    rt_extension: str
    forward_primer: str
    forward_3p_pos: int              # 0-based miRNA position of fwd primer 3' end
    universal_reverse: str
    siblings: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    notes: list = field(default_factory=list)


def design_assay(target_name: str, target_seq: str, family: dict,
                 forward_tail: str = "", min_forward_len: int = 15) -> StemLoopAssay:
    t = to_dna(target_seq)
    n = len(t)
    rt_window_start = n - RT_EXTENSION_NT      # first position inside the RT 6-nt window

    # --- RT primer ---
    rt_extension = revcomp(t[-RT_EXTENSION_NT:])
    rt_primer = UNIVERSAL_BACKBONE + rt_extension

    # --- forward primer ---
    # canonical: covers miRNA up to (but excluding) the last 6 nt; 3' end at n-7.
    fwd_3p = rt_window_start - 1                # 0-based, last position before RT window
    fwd_core = t[:fwd_3p + 1]
    if len(fwd_core) < min_forward_len:
        # too short -> need a 5' tail to reach usable length/Tm
        pass
    forward_primer = forward_tail + fwd_core

    # --- per-sibling discrimination, regime-aware ---
    sibs = {k: to_dna(v) for k, v in family.items() if k != target_name}
    results = []
    for name, s in sibs.items():
        diffs = [i for i in range(min(n, len(s))) if t[i] != s[i]]
        if not diffs:
            results.append(AssaySibling(name, None, None, "identical", None, 99.0,
                                        "Identical over compared region; cross-reacts fully."))
            continue
        nearest = max(diffs)
        from_3p = (n - 1) - nearest

        if nearest <= fwd_3p:
            # reachable by forward-primer 3' end: strong kinetic discrimination,
            # IF we place the forward 3' exactly on it. Report that placement.
            sib_base = s[nearest]
            ev = emp.evaluate(detected=t[nearest], other=sib_base, orientation="forward",
                              deliberate_minus2=(t[nearest - 1] if nearest > 0 else "A"),
                              genomic_minus2_sense=(t[nearest - 1] if nearest > 0 else "A"),
                              rt_rna_template=False)  # qPCR DNA template
            results.append(AssaySibling(
                name, nearest, from_3p, "forward-terminal",
                f"{t[nearest]}.{COMP[sib_base]}", round(ev.discrimination_log, 2),
                f"Place forward primer 3' on position {nearest+1} (qPCR kinetic block)."))
        else:
            # difference falls in the RT 6-nt window: RT-step discrimination only
            sib_base = s[nearest]
            ev = emp.evaluate(detected=t[nearest], other=sib_base, orientation="forward",
                              deliberate_minus2=(t[nearest - 1] if nearest > 0 else "A"),
                              genomic_minus2_sense=(t[nearest - 1] if nearest > 0 else "A"),
                              rt_rna_template=True)   # RNA template, wobble-aware
            # inside the 6-mer it is internal, so even the wobble value is optimistic;
            # treat as weak.
            disc = min(ev.discrimination_log, 1.0)
            results.append(AssaySibling(
                name, nearest, from_3p, "rt-window",
                f"{t[nearest]}.{COMP[sib_base]}", round(disc, 2),
                f"Difference at position {nearest+1} lies in the RT 6-nt anneal window; "
                "discrimination is RT-step only (internal, weak). To exploit a qPCR "
                "kinetic block, shorten the RT extension so the forward primer can reach "
                "this base (specificity/Tm tradeoff)."))

    assay = StemLoopAssay(target_name, t, rt_primer, rt_extension,
                          forward_primer, fwd_3p, UNIVERSAL_REVERSE, siblings=results)

    # warnings
    weak = [r.name for r in results if r.regime in ("rt-window",) and r.discrimination_log < 2.0]
    strong = [r.name for r in results if r.regime == "forward-terminal" and r.discrimination_log >= 2.0]
    ident = [r.name for r in results if r.regime == "identical"]
    if strong:
        assay.notes.append("Cleanly discriminated at the qPCR step (forward-primer 3' "
                           "placement available): " + ", ".join(strong) + ".")
    if weak:
        assay.warnings.append("Weak discrimination (difference sits in the RT 6-nt window; "
                              "RT-step only): " + ", ".join(weak) + ". Residual cross-"
                              "reactivity expected. Options: shorten RT extension to reach "
                              "the base at the forward-primer 3' end, or use an LNA base.")
    if ident:
        assay.warnings.append("Cannot be distinguished at all (identical over compared "
                              "region): " + ", ".join(ident) + ".")
    if tm_wallace(assay.forward_primer) < 50:
        assay.notes.append(f"Forward primer Tm is low (~{tm_wallace(assay.forward_primer)}C "
                           "by Wallace); add a 5' tail (forward_tail=) to raise it.")
    return assay


def assay_report(assay: StemLoopAssay) -> str:
    lines = []
    lines.append(f"=== Stem-loop RT-qPCR assay for {assay.target} ===")
    lines.append(f"miRNA (cDNA sense): 5'-{assay.target_dna}-3'  ({len(assay.target_dna)} nt)")
    lines.append("")
    lines.append("PRIMERS")
    lines.append(f"  RT stem-loop primer : 5'-{assay.rt_primer}-3'")
    lines.append(f"      (universal backbone + 6-nt extension '{assay.rt_extension}' = "
                 f"revcomp of miRNA 3' end)")
    lines.append(f"  Forward primer      : 5'-{assay.forward_primer}-3'  "
                 f"(3' end at miRNA pos {assay.forward_3p_pos+1})")
    lines.append(f"  Universal reverse   : 5'-{assay.universal_reverse}-3'")
    lines.append("")
    lines.append("SIBLING DISCRIMINATION")
    for s in assay.siblings:
        pos = (s.diff_pos + 1) if s.diff_pos is not None else "--"
        tag = {"forward-terminal": "STRONG", "rt-window": "weak ",
               "internal": "none ", "identical": "NONE "}.get(s.regime, "?")
        lines.append(f"  [{tag}] {s.name}: diff@pos {pos} ({s.from_3p} from 3'), "
                     f"mismatch {s.terminal_mismatch}, ~{s.discrimination_log} log  [{s.regime}]")
    if assay.notes:
        lines.append("")
        for nnote in assay.notes:
            lines.append("  note: " + nnote)
    if assay.warnings:
        lines.append("")
        for w in assay.warnings:
            lines.append("  ! " + w)
    return "\n".join(lines)
