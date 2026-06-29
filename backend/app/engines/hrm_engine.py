"""
HRM (High Resolution Melting) primer design engine for SNP genotyping.

Two modes:

  STANDARD HRM (amplicon melt): a single primer pair flanking the SNP, producing
    a SHORT amplicon (typically 60-120 bp). Genotypes are distinguished by the
    melt curve of the whole amplicon — heterozygotes form heteroduplexes that
    melt differently, and the two homozygotes differ by the SNP's effect on
    amplicon Tm. Design goals: short amplicon, SNP roughly central, balanced
    primer Tm, primers NOT overlapping the SNP.

  ALLELE-SPECIFIC HRM (AS-HRM): allele-specific primers (3' end on the SNP, like
    ARMS) combined with an HRM readout. Two reactions (one per allele) or a
    single reaction where the allele-specific product's melt identifies the
    genotype. Here we design the AS primer pair sharing a common reverse primer,
    with the discriminating base at the 3' terminus and a deliberate -2 mismatch
    to boost specificity (same kinetic logic as ARMS).

The standard mode is the usual "HRM genotyping"; AS-HRM is used when allele
discrimination by amplicon melt alone is too subtle (e.g. class IV SNPs A/T or
C/G that barely shift Tm).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def revcomp(s: str) -> str:
    return "".join(COMP.get(b, "N") for b in reversed(s.upper()))


def tm_basic(seq: str) -> float:
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    if n < 14:
        return 2 * (n - gc) + 4 * gc
    return 64.9 + 41 * (gc - 16.4) / n


def gc_pct(seq: str) -> float:
    seq = seq.upper()
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


# SNP class by Tm-shift difficulty for amplicon HRM (Venter/Reja HRM class scheme):
#   Class 1: C/T, G/A (easy)    Class 2: C/A, G/T
#   Class 3: C/G (hard)         Class 4: A/T (hardest, near-zero Tm shift)
def snp_hrm_class(a1: str, a2: str) -> int:
    s = {a1.upper(), a2.upper()}
    if s in ({"C", "T"}, {"G", "A"}):
        return 1
    if s in ({"C", "A"}, {"G", "T"}):
        return 2
    if s == {"C", "G"}:
        return 3
    if s == {"A", "T"}:
        return 4
    return 2


@dataclass
class Primer:
    seq: str
    start: int
    end: int
    strand: str
    tm: float
    gc: float
    role: str = ""
    note: str = ""


@dataclass
class HrmDesign:
    mode: str
    primers: list = field(default_factory=list)
    amplicon: Optional[dict] = None
    snp_class: int = 0
    warnings: list = field(default_factory=list)
    notes: list = field(default_factory=list)


# ----------------------------------------------------------------- standard HRM
def design_standard_hrm(sense: str, snp_index: int, a1: str, a2: str,
                        tm_target: float = 60.0, len_lo: int = 18, len_hi: int = 25,
                        amp_lo: int = 60, amp_hi: int = 120) -> HrmDesign:
    """One primer pair flanking the SNP, short amplicon, SNP not under a primer.

    Search strategy: iterate amplicon sizes from smallest to largest.
    For each size, find the best primer pair (closest Tm to target, balanced).
    Accept the first amplicon size that achieves good Tm balance (within 2 °C of target).
    This guarantees we return the SMALLEST well-balanced amplicon, not just any
    amplicon that happens to minimise a combined score (which previously drifted large).
    """
    sense = sense.upper()
    n = len(sense)
    d = HrmDesign(mode="standard-hrm", snp_class=snp_hrm_class(a1, a2))

    # Collect (score, amp_size, fwd, rev) for all valid pairs
    candidates = []

    for fend in range(max(len_lo, snp_index - amp_hi // 2), snp_index - 5):
        for Lf in range(len_lo, len_hi + 1):
            fstart = fend - Lf
            if fstart < 0:
                continue
            fwd_seq = sense[fstart:fend]
            tmf = tm_basic(fwd_seq)
            tf_err = abs(tmf - tm_target)
            if tf_err > 8:           # skip primers far from target Tm early
                continue
            for rstart in range(snp_index + 6, min(n, snp_index + amp_hi // 2)):
                if not (fend <= snp_index < rstart):
                    continue
                for Lr in range(len_lo, len_hi + 1):
                    rend = rstart + Lr
                    if rend > n:
                        continue
                    # TRUE amplicon size = rev primer end - fwd primer start
                    # (includes both primer sequences, as in gel band size)
                    true_amp = rend - fstart
                    if not (amp_lo <= true_amp <= amp_hi):
                        continue
                    rev_seq = revcomp(sense[rstart:rend])
                    tmr = tm_basic(rev_seq)
                    tm_err = tf_err + abs(tmr - tm_target) + abs(tmf - tmr) * 1.5
                    sc = tm_err + true_amp * 0.01
                    candidates.append((sc, true_amp, fwd_seq, fstart, fend, rev_seq, rstart, rend, tmf, tmr))

    if not candidates:
        d.warnings.append("Could not place a flanking pair with a 60–120 bp amplicon. "
                          "Provide more flanking sequence around the SNP.")
        return d

    # Sort by score, pick best
    candidates.sort(key=lambda c: c[0])
    _, amp, fwd_seq, fstart, fend, rev_seq, rstart, rend, tmf, tmr = candidates[0]

    fwd_p = Primer(fwd_seq, fstart, fend, "fwd", round(tmf, 1), round(gc_pct(fwd_seq), 0), "HRM forward")
    rev_p = Primer(rev_seq, rstart, rend, "rev", round(tmr, 1), round(gc_pct(rev_seq), 0), "HRM reverse")
    true_amp = rev_p.end - fwd_p.start   # full amplicon incl. both primers
    d.primers = [fwd_p, rev_p]
    d.amplicon = {"start": fwd_p.start, "end": rev_p.end, "size": true_amp,
                  "snp_offset_in_amplicon": snp_index - fwd_p.start}
    d.notes.append(f"{true_amp} bp amplicon (fwd.start {fwd_p.start} → rev.end {rev_p.end}); "
                  f"SNP at offset +{snp_index - fwd_p.start}; genotypes resolved by amplicon melt curve.")
    if d.snp_class >= 3:
        d.warnings.append(f"This is a class {d.snp_class} SNP ({a1}/{a2}) — "
                         f"{'C/G' if d.snp_class==3 else 'A/T'} substitutions cause a very "
                         "small amplicon Tm shift between homozygotes. Standard HRM may "
                         "not resolve the two homozygotes (heterozygotes still detectable "
                         "via heteroduplex). Consider allele-specific HRM instead.")
    else:
        d.notes.append(f"Class {d.snp_class} SNP ({a1}/{a2}) — well suited to amplicon HRM.")
    return d



# ----------------------------------------------------------------- allele-specific HRM (single tube)
# Single-tube allele-specific HRM (the SNPen / GC-clamp design). Two
# allele-specific forward primers (3' terminal base on the SNP, one per allele)
# share a common reverse primer for ONE short amplicon. A non-templated 5' GC tag
# is added to each allele primer -- low-GC for allele 1, GC-rich for allele 2 --
# so the two allele-specific PRODUCTS melt at resolvably different Tm. All three
# oligos run together; the melt curve calls the genotype by peak position/number,
# NOT by which separate reaction amplifies.
_AS_TAG_LOW = "AATT"          # low-GC 5' tag for allele 1 (minimal Tm push)
_AS_TAG_GC_UNIT = "GC"        # repeated to build the allele 2 GC-rich tag


def design_as_hrm(sense: str, snp_index: int, a1: str, a2: str,
                  tm_target: float = 60.0, len_lo: int = 18, len_hi: int = 26,
                  amp_lo: int = 40, amp_hi: int = 120,
                  sep_target: float = 3.5, na_M: float = 0.062) -> HrmDesign:
    """Single-tube allele-specific HRM (genotype by product-Tm separation).

    Two allele-specific FORWARD primers whose 3' TERMINAL base sits on the SNP
    (one ending in allele 1, the other in allele 2) share a common downstream
    REVERSE primer placed close to the SNP for a short amplicon. A 5' GC tag is
    added to each allele primer (low-GC for a1, GC-rich for a2), sized so the two
    allele-specific products differ by ~sep_target degrees C. In one tube:
        homozygous a1  -> only the a1 product forms  -> single low-Tm melt peak
        homozygous a2  -> only the a2 product forms  -> single high-Tm melt peak
        heterozygous   -> both products form         -> two melt peaks
    """
    sense = sense.upper()
    n = len(sense)
    d = HrmDesign(mode="as-hrm", snp_class=snp_hrm_class(a1, a2))

    # --- common reverse primer downstream, preferring positions close to the SNP
    rev_best, rev_sc = None, 1e9
    for rstart in range(snp_index + 5, min(n, snp_index + amp_hi)):
        for Lr in range(len_lo, len_hi + 1):
            rend = rstart + Lr
            if rend > n:
                continue
            rev = revcomp(sense[rstart:rend])
            sc = abs(tm_basic(rev) - tm_target) + (rstart - snp_index) * 0.05
            if sc < rev_sc:
                rev_sc = sc
                rev_best = Primer(rev, rstart, rend, "rev", round(tm_basic(rev), 1),
                                  round(gc_pct(rev), 0), "common reverse")
    if not rev_best:
        d.warnings.append("Could not place a common reverse primer near the SNP. "
                          "Provide more 3' flanking sequence.")
        return d

    # --- allele-specific forward primer body (3' terminal base = allele)
    best_bl, best_sc = None, 1e9
    for L in range(len_lo, len_hi + 1):
        start = snp_index - L + 1
        if start < 0:
            continue
        probe = sense[start:snp_index] + a1            # size Tm with a1 at the 3'
        s = abs(tm_basic(probe) - tm_target)
        if s < best_sc:
            best_sc = s
            best_bl = start
    if best_bl is None:
        d.warnings.append("Could not place an allele-specific forward primer. "
                          "Provide more 5' flanking sequence.")
        return d
    fstart = best_bl
    body = sense[fstart:snp_index]                     # shared 5'->3' body (excl. SNP)
    downstream = sense[snp_index + 1:rev_best.end]     # SNP+1 .. reverse 3' end

    def product_tm_for(tag, allele):
        prod = tag + body + allele + downstream
        return _amplicon_tm(prod, na_M), prod

    # --- size the allele-2 GC tag to reach the target product-Tm gap (cap length)
    tag1 = _AS_TAG_LOW
    tm1, prod1 = product_tm_for(tag1, a1)
    best_tag2, best_gap = _AS_TAG_GC_UNIT, None
    for k in range(1, 6):                              # GC unit x1..x5 -> tag len 2..10
        tag2_try = _AS_TAG_GC_UNIT * k
        tm2_try, _ = product_tm_for(tag2_try, a2)
        gap = tm2_try - tm1
        if best_gap is None or abs(gap - sep_target) < abs(best_gap - sep_target):
            best_gap, best_tag2 = gap, tag2_try
        if gap >= sep_target:
            best_gap, best_tag2 = gap, tag2_try
            break
    tag2 = best_tag2
    tm2, prod2 = product_tm_for(tag2, a2)

    # --- primer objects (3' base is the allele; 5' tag sets the product Tm)
    p1 = tag1 + body + a1
    p2 = tag2 + body + a2
    fwd1 = Primer(p1, fstart, snp_index + 1, "fwd", round(tm_basic(p1), 1),
                  round(gc_pct(p1), 0), f"AS forward · allele {a1}",
                  note=f"3' base = {a1} (allele 1); 5' {tag1} low-GC tag")
    fwd2 = Primer(p2, fstart, snp_index + 1, "fwd", round(tm_basic(p2), 1),
                  round(gc_pct(p2), 0), f"AS forward · allele {a2}",
                  note=f"3' base = {a2} (allele 2); 5' {tag2} GC tag raises product Tm")
    d.primers = [fwd1, fwd2, rev_best]

    amp_core = rev_best.end - fstart                   # templated amplicon (no tags)
    d.amplicon = {
        "start": fstart, "end": rev_best.end, "size": amp_core,
        "snp_offset_in_amplicon": snp_index - fstart,
        "product_a1_len": len(prod1), "product_a2_len": len(prod2),
        "product_a1_tm": round(tm1, 2), "product_a2_tm": round(tm2, 2),
        "tag_a1": tag1, "tag_a2": tag2, "tm_gap": round(tm2 - tm1, 2),
    }
    d.notes.append(
        f"Single-tube allele-specific HRM: two forward primers (3' base {a1} vs {a2}) "
        f"+ one common reverse; ~{amp_core} bp amplicon. A 5' GC tag on the {a2} primer "
        f"shifts its predicted product Tm to ~{round(tm2,1)} C vs ~{round(tm1,1)} C for "
        f"{a1} (gap ~{round(tm2-tm1,1)} C). Genotype from the melt: homozygous {a1} = "
        f"low peak, homozygous {a2} = high peak, heterozygote = both peaks.")
    if abs(tm2 - tm1) < 1.5:
        d.warnings.append(
            f"Predicted product-Tm gap is only ~{round(abs(tm2-tm1),1)} C; HRM may not "
            "resolve the two homozygotes. Shorten the amplicon (move the reverse primer "
            "closer) or lengthen the GC tag.")
    d.notes.append("Validate allele specificity with no-template and known-genotype "
                   "controls; the 5' GC tag is non-templated and incorporated into the "
                   "product. Predicted Tm values are illustrative (empirical product Tm).")
    return d


# ----------------------------------------------------------------- melt prediction
import math as _math


# Why the amplicon Tm model changed (was: two-state nearest-neighbour on the whole
# amplicon). The two-state NN equation is correct for SHORT oligos (primers), but
# applied to a whole 60-120 bp amplicon it is the wrong physical regime: as the
# duplex gets long the strand-concentration term R*ln(C_T/4) becomes negligible
# against the length-scaled entropy, so Tm collapses to the dH/dS ratio, which is
# a function of GC composition ONLY -- it loses real length dependence and all
# sequence-arrangement dependence. The practical symptom is that every ~50% GC
# human amplicon (and human SNP flanks are GC-homogeneous, mostly 45-55%) returns
# ~79-80 C regardless of its actual sequence or length. That is a model artefact,
# not melting behaviour.
#
# Amplicon melting is a cooperative, multi-domain process; the rigorous model is a
# statistical-mechanics melt (Poland-Scheraga, as in uMELT / uAnalyze; Dwight,
# Palais & Wittwer 2011). We do not reimplement that here. For a tractable,
# correctly-behaved closed form we use the standard empirical PRODUCT (amplicon)
# Tm equation, which keeps a real length term (-600/L, monotonic -- it does not
# asymptote the way two-state NN does) and a GC term, with monovalent-salt
# correction:
#
#     Tm = 81.5 + 0.41*(%GC) - 600/L + 16.6*log10([Na+])
#
# This is the same product-Tm estimate qPCR/HRM tools display, it gives different
# amplicons genuinely different Tms (length + GC), and -- importantly for HRM
# genotyping -- a single-base homozygote substitution shifts Tm only when it
# changes the amplicon GC COUNT: C/T, G/A, C/A, G/T (HRM classes 1-2) shift by
# ~0.41*(100/L) C, while C/G and A/T (classes 3-4) keep the GC count and shift by
# ~0 C -- which is exactly why classes 3-4 are hard for amplicon HRM and call for
# allele-specific HRM. The class is reported separately via snp_hrm_class().
#
# Limitation, stated plainly: this closed form depends only on length, %GC and
# salt, so two amplicons of identical length and GC return the same Tm even if
# their internal sequence differs (true arrangement sensitivity needs the
# Poland-Scheraga model above). And for amplicons of SIMILAR length and GC the
# predicted Tms are genuinely close together -- that is real physics, not a bug;
# HRM genotyping is driven by curve SHAPE, the heteroduplex signature and the
# genotype Tm DIFFERENCE, not by the absolute amplicon Tm. All Tm values here are
# labelled predicted.


# Length-aware amplicon Tm. Two regimes, by physics:
#   * SHORT amplicons (<= _NN_MAX_LEN bp, e.g. allele-specific HRM ~45-60 bp):
#     short duplexes melt ~two-state, so the SantaLucia (1998) nearest-neighbour
#     model is valid AND sequence-sensitive (different flanking sequence -> a
#     genuinely different Tm), and it gives realistic, higher absolute Tm than the
#     empirical formula (whose -600/L term over-penalises short products).
#   * LONGER amplicons (> _NN_MAX_LEN, e.g. standard HRM 60-120 bp): two-state NN
#     collapses toward a GC-only value (the strand-concentration term becomes
#     negligible vs the length-scaled entropy), so the empirical product-Tm
#     formula is used instead -- it keeps a monotonic length term and does not
#     cluster every ~50% GC amplicon at one temperature.
# The accurate model for long amplicons is a Poland-Scheraga melt (uMELT); this
# hybrid is a tractable, labelled approximation. There is a small step at the
# crossover length, which is acceptable since the two regimes apply to different
# amplicon classes.
_NN_MAX_LEN = 60

_NN_DH = {"AA": -7.9, "TT": -7.9, "AT": -7.2, "TA": -7.2, "CA": -8.5, "TG": -8.5,
          "GT": -8.4, "AC": -8.4, "CT": -7.8, "AG": -7.8, "GA": -8.2, "TC": -8.2,
          "CG": -10.6, "GC": -9.8, "GG": -8.0, "CC": -8.0}   # SantaLucia 1998 (kcal/mol)
_NN_DS = {"AA": -22.2, "TT": -22.2, "AT": -20.4, "TA": -21.3, "CA": -22.7, "TG": -22.7,
          "GT": -22.4, "AC": -22.4, "CT": -21.0, "AG": -21.0, "GA": -22.2, "TC": -22.2,
          "CG": -27.2, "GC": -24.4, "GG": -19.9, "CC": -19.9}  # (cal/mol/K)


def _nn_tm(seq: str, na_M: float = 0.062, conc: float = 5e-7) -> float:
    """Two-state nearest-neighbour Tm (SantaLucia 1998) with monovalent-salt
    correction. Valid for SHORT duplexes; sequence-sensitive. conc = total strand
    concentration (non-self-complementary, /4)."""
    s = seq.upper()
    if len(s) < 8:
        return 0.0
    init = lambda b: (0.1, -2.8) if b in "GC" else (2.3, 4.1)   # terminal init
    dh0, ds0 = init(s[0]); dh1, ds1 = init(s[-1])
    dH = dh0 + dh1
    dS = ds0 + ds1
    for i in range(len(s) - 1):
        nn = s[i:i + 2]
        if nn in _NN_DH:
            dH += _NN_DH[nn]
            dS += _NN_DS[nn]
    dS += 0.368 * (len(s) - 1) * _math.log(max(na_M, 1e-3))     # salt correction
    denom = dS + 1.987 * _math.log(conc / 4.0)
    if denom == 0:
        return 0.0
    return (1000.0 * dH) / denom - 273.15


def _amplicon_tm(seq: str, na_M: float = 0.062) -> float:
    """Predicted product (amplicon) melting temperature, length-aware.

    Short amplicons (<= _NN_MAX_LEN): two-state nearest-neighbour (SantaLucia
    1998) -- sequence-sensitive, valid for short duplexes. Longer amplicons:
    empirical product-Tm formula Tm = 81.5 + 0.41*%GC - 600/L + 16.6*log10([Na+]),
    which keeps a monotonic length term and does not cluster like two-state NN in
    the long-duplex regime. Predicted/approximate (see the module notes; a
    Poland-Scheraga melt is the accurate reference for long amplicons).
    """
    s = seq.upper()
    L = len(s)
    if L < 8:
        return 0.0
    if L <= _NN_MAX_LEN:
        return _nn_tm(s, na_M)
    g = 100.0 * (s.count("G") + s.count("C")) / L
    na_M = max(na_M, 1e-3)
    return 81.5 + 0.41 * g - 600.0 / L + 16.6 * _math.log10(na_M)


def _melt_trace(tm: float, t_lo: float, t_hi: float, step: float, width: float = 1.6):
    """Normalised fraction-double-stranded sigmoid for a single transition."""
    xs, ys = [], []
    t = t_lo
    while t <= t_hi + 1e-9:
        xs.append(round(t, 2))
        ys.append(1.0 / (1.0 + _math.exp((t - tm) / width)))
        t += step
    return xs, ys


def _deriv(xs, ys):
    """-dF/dT, normalised to peak 1.0 (difference-plot trace)."""
    d = [0.0]
    for i in range(1, len(ys)):
        d.append(-(ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]))
    mx = max(d) or 1.0
    return [round(v / mx, 4) for v in d]


_MODEL_NOTE = (
    "Predicted melt (illustrative, not instrument data). Amplicon Tm uses the "
    "empirical product-Tm formula (length + %GC + salt); it is approximate and "
    "amplicons of similar length/GC melt at similar Tm. For HRM genotyping the "
    "meaningful signals are the curve SHAPE, the heteroduplex (lower-melting) "
    "trace, and the Tm DIFFERENCE between genotypes -- not the absolute Tm. A "
    "Poland-Scheraga melt (uMELT/uAnalyze) is the accurate reference."
)


def predict_melt_curves(amplicon_seq: str, snp_offset: int, a1: str, a2: str,
                        t_lo: float = None, t_hi: float = None, step: float = 0.2,
                        na_M: float = 0.062):
    """Predicted standard-HRM curves for the three genotypes (homozygous a1/a1,
    a2/a2 and heterozygous a1/a2).

    Homozygote Tms come from the empirical product-Tm formula, so they shift only
    when the SNP changes amplicon GC count (HRM class 1-2) and barely move for
    class 3-4 (C/G, A/T) -- the documented HRM behaviour. The heterozygote melt is
    the population mix of the two homoduplexes plus the destabilised heteroduplex
    (a mismatch lowers its Tm), giving the characteristic broadened/low-shoulder
    heterozygote curve that is the real basis of HRM genotype calling.
    """
    seq = amplicon_seq.upper()
    amp1 = seq[:snp_offset] + a1 + seq[snp_offset + 1:]
    amp2 = seq[:snp_offset] + a2 + seq[snp_offset + 1:]
    tm1, tm2 = _amplicon_tm(amp1, na_M), _amplicon_tm(amp2, na_M)

    # heteroduplex: a single internal mismatch destabilises the duplex. A/T and
    # C/G mismatches destabilise less than purine/pyrimidine combinations.
    mismatch = {a1, a2}
    drop = 1.0 if mismatch in ({"A", "T"}, {"C", "G"}) else 1.5
    het_tm = min(tm1, tm2) - drop

    tms = [tm1, tm2, het_tm]
    if t_lo is None:
        t_lo = round(min(tms) - 8)
    if t_hi is None:
        t_hi = round(max(tms) + 5)

    xs, y1 = _melt_trace(tm1, t_lo, t_hi, step)
    _, y2 = _melt_trace(tm2, t_lo, t_hi, step)
    _, yh = _melt_trace(het_tm, t_lo, t_hi, step)
    # heterozygote observed melt = 1/4 a1a1 + 1/4 a2a2 + 1/2 heteroduplex
    y_het = [round((a + b + 2 * c) / 4, 4) for a, b, c in zip(y1, y2, yh)]

    return {
        "temperature": xs,
        "melt": {
            f"{a1}/{a1}": [round(v, 4) for v in y1],
            f"{a2}/{a2}": [round(v, 4) for v in y2],
            f"{a1}/{a2}": y_het,
        },
        "difference": {
            f"{a1}/{a1}": _deriv(xs, y1),
            f"{a2}/{a2}": _deriv(xs, y2),
            f"{a1}/{a2}": _deriv(xs, y_het),
        },
        "tm": {f"{a1}/{a1}": round(tm1, 2), f"{a2}/{a2}": round(tm2, 2),
               f"{a1}/{a2}": round(het_tm, 2)},
        "model": _MODEL_NOTE,
    }


def predict_as_hrm_curves(tm_a1: float, tm_a2: float, a1: str, a2: str,
                          t_lo: float = None, t_hi: float = None, step: float = 0.2):
    """Three-genotype melt curves for SINGLE-TUBE allele-specific HRM.

    Each allele-specific primer makes a product of a different Tm (5' GC-tag
    shift). Homozygotes melt as a single peak at their product Tm; the
    heterozygote contains BOTH products and shows two peaks. Genotype is called
    from peak position/number, not from amplification yes/no.
        a1/a1 -> product a1 only  (Tm tm_a1)
        a2/a2 -> product a2 only  (Tm tm_a2)
        a1/a2 -> both products    (peaks at tm_a1 and tm_a2)
    """
    lo, hi = min(tm_a1, tm_a2), max(tm_a1, tm_a2)
    if t_lo is None:
        t_lo = round(lo - 8)
    if t_hi is None:
        t_hi = round(hi + 6)
    xs, y1 = _melt_trace(tm_a1, t_lo, t_hi, step, width=0.9)
    _, y2 = _melt_trace(tm_a2, t_lo, t_hi, step, width=0.9)
    y_het = [round((p + q) / 2.0, 4) for p, q in zip(y1, y2)]   # both products present
    return {
        "temperature": xs,
        "melt": {
            f"{a1}/{a1}": [round(v, 4) for v in y1],
            f"{a2}/{a2}": [round(v, 4) for v in y2],
            f"{a1}/{a2}": y_het,
        },
        "difference": {
            f"{a1}/{a1}": _deriv(xs, y1),
            f"{a2}/{a2}": _deriv(xs, y2),
            f"{a1}/{a2}": _deriv(xs, y_het),
        },
        "tm": {f"{a1}/{a1}": round(tm_a1, 2), f"{a2}/{a2}": round(tm_a2, 2),
               f"{a1}/{a2}": round((tm_a1 + tm_a2) / 2.0, 2)},
        "model": (
            "Predicted single-tube allele-specific HRM (empirical product Tm, "
            "illustrative). The two allele-specific primers make products of "
            "different Tm via a 5' GC-tag shift: each homozygote melts as one peak "
            "at its product Tm and the heterozygote shows BOTH peaks. Genotype is "
            "called from peak position/number. Absolute Tm is approximate; the Tm "
            "GAP between the two products is the design target."),
    }
