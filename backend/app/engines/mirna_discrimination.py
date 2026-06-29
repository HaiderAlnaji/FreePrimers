"""
miRNA cross-family discrimination engine.

Problem: miRNA families (e.g. let-7) contain members differing by as little as
one nucleotide. A qPCR assay for one member cross-reacts with its siblings.
Published stem-loop / TaqMan assays achieve only 0.1-3.7% residual
cross-reactivity for single-base differences, and the worst offenders are
G.U (G.T) wobble differences left INTERNAL to the RT duplex (Chen et al.
2005, NAR: "most cross-reactions resulted from G-T mismatches during the RT
reaction, let-7a vs let-7c").

Key idea (transfers problem #1's kinetic model):
  An internal difference discriminates weakly (wobble, near-zero kinetic
  block). But if the forward qPCR primer's 3' TERMINUS is placed ON the
  discriminating base, the difference becomes a 3'-terminal mismatch during
  qPCR -- which the polymerase resists extending (Huang 1992 extension-rate
  model). A deliberate -2 mismatch can be stacked on top (ARMS trick).

  When the discriminating base is too far from the miRNA 3' end to sit at a
  primer terminus (e.g. let-7f differs at position 12, mid-sequence), NO
  primer-3' placement can reach it -> the engine reports it as
  undistinguishable by primer design, consistent with the patent literature
  ("no amplification selection during qPCR" for internal mismatches).

This module:
  - aligns family members to the target (equal-length mature miRNAs),
  - finds every position where a sibling differs,
  - for each sibling, finds the best forward-primer 3' placement that turns a
    near-3' difference into a terminal (or near-terminal) mismatch,
  - scores residual cross-reactivity with the empirical kinetic model,
  - flags undistinguishable siblings honestly.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from . import empirical_discrimination as emp

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def to_dna(seq: str) -> str:
    return seq.strip().upper().replace("U", "T")


@dataclass
class SiblingResult:
    name: str
    diff_positions: list          # 0-based positions where sibling != target (aligned)
    nearest_3p_diff: Optional[int]  # position closest to 3' end
    from_3p: Optional[int]         # how many nt from the 3' end that diff is
    best_primer_3p: Optional[int]  # chosen forward-primer 3' position (0-based)
    terminal_mismatch: Optional[str]  # e.g. "A.C" primer.template at terminus
    discrimination_log: float      # predicted log-fold suppression (higher=better)
    distinguishable: bool
    note: str = ""
    required_primer_3p: Optional[int] = None  # primer 3' pos this sibling needs


@dataclass
class MirnaDesign:
    target: str
    target_seq: str
    primer_3p_pos: int             # 0-based position of forward-primer 3' end
    primer_len: int
    forward_primer: str
    deliberate_minus2: Optional[str]
    siblings: list = field(default_factory=list)
    warnings: list = field(default_factory=list)


def align_family(target_seq: str, family: dict) -> dict:
    """Return {name: [diff positions]} for each sibling vs target.

    Assumes equal-length mature miRNAs (true within most miRNA families after
    5' alignment). Unequal lengths are handled by comparing over the shared
    5'-anchored length and flagging the overhang.
    """
    t = to_dna(target_seq)
    out = {}
    for name, seq in family.items():
        s = to_dna(seq)
        n = min(len(t), len(s))
        diffs = [i for i in range(n) if t[i] != s[i]]
        if len(t) != len(s):
            diffs.append(("len", len(t), len(s)))
        out[name] = diffs
    return out


def _terminal_mismatch_for_sibling(target_dna, sibling_dna, primer_3p_pos):
    """At a given forward-primer 3' position, what mismatch (if any) does the
    sibling present at the primer terminus?

    Forward primer is the target sense up to primer_3p_pos (inclusive). In qPCR
    it anneals to the antisense cDNA. For the sibling template, the base at
    primer_3p_pos on the antisense strand is COMP[sibling_sense[pos]]. The
    primer 3' base is target_sense[pos]. Mismatch iff
    primer_base != COMP[antisense_base] = sibling_base.
    Returns (primer_base, sibling_template_base_on_antisense, is_mismatch).
    """
    primer_base = target_dna[primer_3p_pos]
    sib_sense = sibling_dna[primer_3p_pos] if primer_3p_pos < len(sibling_dna) else None
    if sib_sense is None:
        return primer_base, None, True
    antisense_base = COMP[sib_sense]
    is_mm = (COMP[primer_base] != antisense_base)  # primer pairs with antisense
    return primer_base, antisense_base, is_mm


def design(target_name: str, target_seq: str, family: dict,
           primer_len_range=(18, 24), allow_deliberate_minus2=True) -> MirnaDesign:
    """Design forward qPCR primer placement(s) to discriminate `target` from
    each sibling.

    A single primer 3' position can only place ONE position at its terminus, so
    siblings whose nearest difference sits at different positions need different
    primers. This function:
      - for each sibling, finds the best reachable primer 3' placement that puts
        that sibling's nearest difference at (or nearest to) the terminus,
      - scores discrimination there with the kinetic model (RT-wobble aware),
      - picks, as the PRIMARY primer, the 3' position that maximises the number
        of siblings discriminated and the minimum discrimination among them,
      - reports per-sibling which primer 3' position each one needs, and flags
        siblings that NO placement can reach (internal differences).
    """
    t = to_dna(target_seq)
    n = len(t)

    sibs = {k: to_dna(v) for k, v in family.items() if k != target_name}
    # nearest-3' single-base difference position for each sibling
    sib_diffpos = {}
    for name, s in sibs.items():
        diffs = [i for i in range(min(len(t), len(s))) if t[i] != s[i]]
        sib_diffpos[name] = max(diffs) if diffs else None

    # candidate primer 3' positions = every sibling's nearest diff that is
    # reachable by a primer within the length range
    cand_positions = sorted({p for p in sib_diffpos.values()
                             if p is not None and primer_len_range[0] - 1 <= p <= n - 1},
                            reverse=True)

    def score_at(p3):
        """Return (sib_results, n_ok, min_disc) for a primer ending at p3."""
        results = _score_all_siblings(target_name, t, sibs, p3, allow_deliberate_minus2)
        ok = [r for r in results if r.distinguishable and r.discrimination_log >= 2.0]
        min_disc = min((r.discrimination_log for r in ok), default=0.0)
        return results, len(ok), min_disc

    # pick PRIMARY placement: maximise (# discriminated, then min discrimination)
    best = None
    for p3 in cand_positions:
        results, n_ok, min_disc = score_at(p3)
        key = (n_ok, min_disc)
        if best is None or key > best[0]:
            best = (key, p3, results)
    if best is None:
        return MirnaDesign(target_name, t, n - 1, n, t, None, siblings=[],
                           warnings=["No reachable single-base differences in family."])

    _, p3, results = best
    L = min(primer_len_range[1], p3 + 1)
    L = max(primer_len_range[0], L)
    start = max(0, p3 - L + 1)
    primer = t[start:p3 + 1]
    minus2 = _pick_minus2(target_name, t, sibs, p3) if allow_deliberate_minus2 and p3 > 0 else None

    # For each sibling, also compute the primer 3' position IT would need (so the
    # report can recommend separate assays where needed).
    for r in results:
        dp = sib_diffpos.get(r.name)
        r.required_primer_3p = dp  # may differ from the primary p3

    warnings = []
    # Classify each sibling by its OWN best achievable outcome, not the primary
    # primer's coverage. A sibling is: discriminated (some reachable primer puts
    # its diff at a terminus with disc>=2), wobble-limited (terminal but G.U),
    # or unreachable (diff too internal for any terminus).
    needs_other, wobble_limited, unreachable = [], [], []
    for r in results:
        dp = r.required_primer_3p
        if dp is None:
            continue
        from_3p_dp = (len(t) - 1) - dp
        # can a primer end exactly on dp?  yes if within length range and dp<n
        reachable = (primer_len_range[0] - 1 <= dp <= len(t) - 1)
        # best-case discrimination if primer ended on dp
        sib_seq = sibs[r.name]
        if dp < len(sib_seq):
            ev = emp.evaluate(detected=t[dp], other=sib_seq[dp], orientation="forward",
                              deliberate_minus2=(t[dp - 1] if dp > 0 else "A"),
                              genomic_minus2_sense=(t[dp - 1] if dp > 0 else "A"),
                              rt_rna_template=True)
            best_disc = ev.discrimination_log
        else:
            best_disc = 0.0
        if not reachable:
            unreachable.append(r.name)
        elif best_disc < 2.0:
            wobble_limited.append((r.name, round(best_disc, 1)))
        elif dp != p3:
            needs_other.append((r.name, dp + 1))

    if needs_other:
        warnings.append("One primer cannot place every sibling's difference at its "
                        "terminus. Separate primers needed: " +
                        "; ".join(f"{nm} (3' on position {pos})" for nm, pos in needs_other) + ".")
    if wobble_limited:
        warnings.append("Weak discrimination even at terminal placement (G.U wobble or "
                        "easy-extension mispair): " +
                        "; ".join(f"{nm} (~{d} log)" for nm, d in wobble_limited) +
                        ". Residual cross-reactivity expected; consider LNA or a distinct "
                        "target region.")
    if unreachable:
        warnings.append("Undistinguishable by any primer 3' placement (difference too far "
                        "from the 3' end; no qPCR amplification selection): " +
                        ", ".join(unreachable) + ".")

    return MirnaDesign(target_name, t, p3, len(primer), primer, minus2,
                       siblings=results, warnings=warnings)


def _score_all_siblings(target_name, t_dna, family, primer_3p_pos, allow_minus2):
    results = []
    for name, seq in family.items():
        if name == target_name:
            continue
        s = to_dna(seq)
        diffs = [i for i in range(min(len(t_dna), len(s))) if t_dna[i] != s[i]]
        nearest = max(diffs) if diffs else None
        if nearest is None:
            results.append(SiblingResult(name, [], None, None, None, None, 99.0, True,
                                         "Identical over compared region."))
            continue
        from_3p = (len(t_dna) - 1) - nearest

        # Is the nearest difference at or within reach of this primer terminus?
        # Discrimination via qPCR requires the diff at/very near the 3' terminus.
        pbase, atemplate, is_mm = _terminal_mismatch_for_sibling(t_dna, s, primer_3p_pos)
        diff_at_terminus = (nearest == primer_3p_pos)

        if diff_at_terminus and is_mm:
            # terminal mismatch -> kinetic block. primer base vs sibling template.
            # template base presented on antisense = COMP[sibling_sense]; mismatch
            # type for the empirical model is primer_base . sibling_sense_base
            sib_base = s[primer_3p_pos]
            ev = emp.evaluate(detected=t_dna[primer_3p_pos], other=sib_base,
                              orientation="forward",
                              deliberate_minus2=(t_dna[primer_3p_pos - 1] if primer_3p_pos > 0 else "A"),
                              genomic_minus2_sense=(t_dna[primer_3p_pos - 1] if primer_3p_pos > 0 else "A"),
                              rt_rna_template=True)
            disc = ev.discrimination_log
            results.append(SiblingResult(
                name, diffs, nearest, from_3p, primer_3p_pos,
                f"{t_dna[primer_3p_pos]}.{COMP[sib_base]}", round(disc, 2), True,
                "Difference placed at primer 3' terminus (kinetic block)."))
        elif from_3p <= 2:
            # near-3' but not at terminus: weak (penultimate) discrimination
            results.append(SiblingResult(
                name, diffs, nearest, from_3p, primer_3p_pos, None, 1.0, True,
                f"Difference {from_3p} nt from 3' end but not at this primer's terminus; "
                "weak (penultimate-zone) discrimination."))
        else:
            # internal difference -> no qPCR selection, only weak RT-step wobble
            results.append(SiblingResult(
                name, diffs, nearest, from_3p, None, None, 0.0, False,
                f"Difference is {from_3p} nt from 3' end (internal); no primer 3' "
                "placement can reach it -> undistinguishable by qPCR selection."))
    return results


def _pick_minus2(target_name, t_dna, family, primer_3p_pos):
    """Pick a deliberate -2 base that raises minimum sibling discrimination."""
    template_m2 = t_dna[primer_3p_pos - 1]
    best = (None, -1.0)
    for d2 in "ACGT":
        if d2 == template_m2:
            continue  # must mismatch target template minimally? keep matched allowed
        # evaluate min discrimination across siblings that sit at terminus
        mins = []
        for name, seq in family.items():
            if name == target_name:
                continue
            s = to_dna(seq)
            diffs = [i for i in range(min(len(t_dna), len(s))) if t_dna[i] != s[i]]
            if not diffs or max(diffs) != primer_3p_pos:
                continue
            ev = emp.evaluate(detected=t_dna[primer_3p_pos], other=s[primer_3p_pos],
                              orientation="forward", deliberate_minus2=d2,
                              genomic_minus2_sense=template_m2)
            mins.append(ev.discrimination_log)
        if mins and min(mins) > best[1]:
            best = (d2, min(mins))
    return best[0]
