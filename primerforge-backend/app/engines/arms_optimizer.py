"""
Joint tetra-primer ARMS-PCR optimizer.

Combines the three validated components on the axes each is valid for:
  - allele discrimination  -> empirical extension-rate model (Huang/Simsek)
  - primer quality          -> primer3 nearest-neighbor thermodynamics
  - band resolvability      -> explicit amplicon-size geometry

Tetra-primer ARMS geometry (Ye et al. 2001), sense coordinates, SNP at S:
  FO  forward outer  : 5' at a  (a < S)            -> control & allele-2 amplicons
  FI  forward inner  : 3' at S, detects one allele -> pairs with RO (downstream allele amplicon)
  RI  reverse inner  : 3' at S, detects other allele -> pairs with FO (upstream allele amplicon)
  RO  reverse outer  : 5' at b  (b > S)            -> control & allele-1 amplicons

Amplicon lengths (= reverse-primer 5' pos - forward-primer 5' pos + 1):
  control  (FO+RO) = b - a + 1
  allele-2 (FO+RI) = (S - a) + Lri        # upstream side; RI detects allele2
  allele-1 (FI+RO) = (b - S) + Lfi        # downstream side; FI detects allele1
Asymmetric outer placement (S-a != b-S) makes the two allele bands resolvable.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import primer3

from . import empirical_discrimination as emp

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
COND = dict(mv_conc=50.0, dv_conc=3.0, dntp_conc=0.8, dna_conc=250.0)


def revcomp(s: str) -> str:
    return "".join(COMP[b] for b in reversed(s))


def gc(s: str) -> float:
    return 100.0 * sum(c in "GC" for c in s) / len(s) if s else 0.0


def tm(s: str) -> float:
    return primer3.calc_tm(s, **COND)


def hairpin_dg(s: str) -> float:
    r = primer3.calc_hairpin(s, **COND)
    return r.dg / 1000.0 if r.structure_found else 0.0


def heterodimer_dg(a: str, b: str) -> float:
    r = primer3.calc_heterodimer(a, b, **COND)
    return r.dg / 1000.0 if r.structure_found else 0.0


# ---- design parameters -----------------------------------------------------
@dataclass
class Params:
    tm_target: float = 60.0
    tm_tol: float = 3.0          # acceptable spread around target
    len_min: int = 18
    len_max: int = 27
    gc_min: float = 30.0
    gc_max: float = 70.0
    allele_size_min: int = 90    # smallest allele amplicon
    allele_size_max: int = 320   # largest allele amplicon
    min_band_sep: int = 40       # min bp between the two allele bands
    min_control_gap: int = 25    # control must exceed larger allele band by this
    dimer_floor: float = -9.0    # cross-dimer dG more negative than this is penalised


# ---- weights for the joint objective (documented; tune via benchmark) ------
W = dict(discrimination=1.0, tm_balance=1.4, gc=0.5, hairpin=0.7,
         dimer=1.0, band=1.2, primability=1.5, size_range=0.8)


@dataclass
class Primer:
    name: str
    seq: str
    role: str                    # 'inner'/'outer'
    tm: float = 0.0
    gc: float = 0.0
    hairpin: float = 0.0
    detail: dict = field(default_factory=dict)


@dataclass
class Design:
    primers: dict                # name -> Primer
    sizes: dict                  # 'control'/'allele1'/'allele2' -> int
    discrimination: dict         # 'FI'/'RI' -> log fold suppression
    score: float = 0.0
    penalties: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)


def _inner_forward(sense: str, S: int, detected: str, other: str, Lfi: int):
    """Build forward inner primer (3' on SNP) with optimised deliberate -2."""
    g2 = sense[S - 1]  # genomic -2 sense base
    best, _ = emp.best_design(detected, other, g2)
    # force forward orientation for FI (it is by construction a forward primer)
    fwd = [c for c in [best] if c.orientation == "forward"]
    if not fwd:
        best, allc = emp.best_design(detected, other, g2)
        fwd = sorted([c for c in allc if c.orientation == "forward"],
                     key=lambda c: (c.discrimination_log, -c.matched_loss_log), reverse=True)
    pick = fwd[0]
    body = sense[S - Lfi + 1:S - 1]       # 5'..-3 (up to -2 exclusive)
    seq = body + pick.deliberate_minus2 + detected
    return seq, pick


def _inner_reverse(sense: str, S: int, detected: str, other: str, Lri: int):
    """Build reverse inner primer (3' on SNP, other strand)."""
    g2 = sense[S + 1]                      # sense base adjacent on the 3' side
    _, allc = emp.best_design(detected, other, g2)
    rev = sorted([c for c in allc if c.orientation == "reverse"],
                 key=lambda c: (c.discrimination_log, -c.matched_loss_log), reverse=True)
    pick = rev[0]
    window = sense[S:S + Lri]              # sense S..S+Lri-1
    rc = revcomp(window)                   # 5'->3'; rc[-1] pairs with sense[S]
    seq = rc[:-2] + pick.deliberate_minus2 + COMP[detected]
    return seq, pick


def _best_outer(sense: str, anchor5: int, forward: bool, p: Params,
                tm_aim: float = None):
    """Pick an outer primer near a 5' anchor, optimising Tm toward a target.

    forward: build sense primer starting at anchor5 (5' end).
    reverse: build reverse primer with 5' end at anchor5 (sense coord).
    tm_aim : Tm to optimise toward. Defaults to p.tm_target, but on GC-locked
             loci the caller passes the (unavoidably high) inner-primer Tm so
             the outers run hot to MATCH the inners rather than sitting at the
             global target and creating a large spread. Also relax the GC
             ceiling when aiming hot, since a high Tm requires high GC.
    Returns (seq, used_5prime_pos, length, tm, gc) or None.
    """
    aim = p.tm_target if tm_aim is None else tm_aim
    gc_ceiling = p.gc_max if aim <= p.tm_target + 1 else min(85.0, p.gc_max + 15)
    best = None
    for L in range(p.len_min, p.len_max + 1):
        for off in range(-4, 5):
            if forward:
                a = anchor5 + off
                if a < 0 or a + L > len(sense):
                    continue
                seq = sense[a:a + L]
                pos5 = a
            else:
                b = anchor5 + off
                if b - L + 1 < 0 or b >= len(sense):
                    continue
                seq = revcomp(sense[b - L + 1:b + 1])
                pos5 = b
            t = tm(seq)
            g = gc(seq)
            if not (p.gc_min <= g <= gc_ceiling):
                continue
            cost = abs(t - aim) + 0.3 * abs(hairpin_dg(seq))
            if best is None or cost < best[0]:
                best = (cost, seq, pos5, L, t, g)
    if best is None:
        return None
    _, seq, pos5, L, t, g = best
    return seq, pos5, L, t, g


def _score(d: Design, primers: dict, sizes: dict, disc: dict, p: Params) -> Design:
    pen = {}
    # Tm balance across all four
    tms = [pr.tm for pr in primers.values()]
    pen["tm_balance"] = (max(tms) - min(tms))
    # individual GC / hairpin
    pen["gc"] = sum(max(0, p.gc_min - pr.gc) + max(0, pr.gc - p.gc_max) for pr in primers.values())
    pen["hairpin"] = sum(max(0.0, -pr.hairpin - 2.0) for pr in primers.values())
    # cross-dimers over all 6 pairs
    names = list(primers)
    worst = 0.0
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            dg = heterodimer_dg(primers[names[i]].seq, primers[names[j]].seq)
            worst = min(worst, dg)
    pen["dimer"] = max(0.0, -(worst - p.dimer_floor))
    # discrimination: want BOTH inners high; penalise the weaker one
    weak = min(disc["FI"], disc["RI"])
    pen["discrimination"] = max(0.0, 6.0 - weak)        # 6 logs ~ ideal
    # primability: matched-allele loss budget already filtered; penalise residual
    pen["primability"] = 0.0
    # band resolvability
    sep = abs(sizes["allele1"] - sizes["allele2"])
    pen["band"] = max(0.0, p.min_band_sep - sep)
    # size range
    inrange = all(p.allele_size_min <= sizes[k] <= p.allele_size_max for k in ("allele1", "allele2"))
    pen["size_range"] = 0.0 if inrange else 30.0
    ctrl_gap = sizes["control"] - max(sizes["allele1"], sizes["allele2"])
    if ctrl_gap < p.min_control_gap:
        pen["size_range"] += (p.min_control_gap - ctrl_gap)

    score = sum(W[k] * v for k, v in pen.items())
    d.penalties = {k: round(v, 2) for k, v in pen.items()}
    d.score = round(score, 2)
    return d


def optimize(sense: str, snp_index: int, allele1: str, allele2: str,
             p: Params = Params(), n_results: int = 3):
    """Design tetra-primer ARMS primers for a SNP.

    sense       : sense-strand sequence (ACGT).
    snp_index   : 0-based index of the SNP base in `sense`.
    allele1/2   : the two alleles (sense strand).
    Returns a list of ranked Design objects (best first).
    """
    sense = sense.upper()
    S = snp_index
    designs = []

    # choose target allele-band size pairs that are resolvable & in range
    size_pairs = []
    lo, hi = p.allele_size_min + 10, p.allele_size_max - 10
    for small in range(lo, hi, 30):
        for big in range(small + p.min_band_sep + 10, hi + 1, 30):
            size_pairs.append((small, big))

    for Lfi in range(p.len_min, p.len_max + 1, 2):
        for Lri in range(p.len_min, p.len_max + 1, 2):
            if S - Lfi + 1 < 0 or S + Lri > len(sense):
                continue
            # search BOTH allele->inner assignments; the better one wins on score.
            for fi_allele, ri_allele in ((allele1, allele2), (allele2, allele1)):
                fi_seq, fi_pick = _inner_forward(sense, S, fi_allele, ri_allele, Lfi)
                ri_seq, ri_pick = _inner_reverse(sense, S, ri_allele, fi_allele, Lri)

                # Inner Tms are fixed by the SNP neighbourhood. If they are
                # forced well above the global target (GC-locked locus), aim the
                # OUTER primers at the inner Tm so all four match, rather than
                # anchoring outers low and creating a large spread.
                inner_tm = (tm(fi_seq) + tm(ri_seq)) / 2.0
                outer_aim = p.tm_target if inner_tm <= p.tm_target + 2 else inner_tm

                for (small, big) in size_pairs:
                    # allele2 (FO+RI) = (S-a)+Lri  ; allele1 (FI+RO) = (b-S)+Lfi
                    # assign the smaller target to one side; try both orientations
                    for sizeA2, sizeB1 in ((small, big), (big, small)):
                        a = S - (sizeA2 - Lri)
                        b = S + (sizeB1 - Lfi)
                        if a < 0 or b >= len(sense):
                            continue
                        fo = _best_outer(sense, a, True, p, outer_aim)
                        ro = _best_outer(sense, b, False, p, outer_aim)
                        if not fo or not ro:
                            continue
                        fo_seq, fo_pos, Lfo, fo_tm, fo_gc = fo
                        ro_seq, ro_pos, Lro, ro_tm, ro_gc = ro

                        sizes = {
                            "control": ro_pos - fo_pos + 1,
                            "allele2": (S - fo_pos) + Lri,
                            "allele1": (ro_pos - S) + Lfi,
                        }
                        primers = {
                            "FO": Primer("FO", fo_seq, "outer", fo_tm, fo_gc, hairpin_dg(fo_seq)),
                            "FI": Primer("FI", fi_seq, "inner", tm(fi_seq), gc(fi_seq), hairpin_dg(fi_seq),
                                         {"detects": fi_allele, "minus2": fi_pick.deliberate_minus2}),
                            "RI": Primer("RI", ri_seq, "inner", tm(ri_seq), gc(ri_seq), hairpin_dg(ri_seq),
                                         {"detects": ri_allele, "minus2": ri_pick.deliberate_minus2}),
                            "RO": Primer("RO", ro_seq, "outer", ro_tm, ro_gc, hairpin_dg(ro_seq)),
                        }
                        disc = {"FI": fi_pick.discrimination_log, "RI": ri_pick.discrimination_log}
                        d = Design(primers, sizes, disc)
                        d.warnings = list(dict.fromkeys(fi_pick.warnings + ri_pick.warnings))
                        _score(d, primers, sizes, disc, p)
                        designs.append(d)

    designs.sort(key=lambda d: d.score)
    return designs[:n_results]
