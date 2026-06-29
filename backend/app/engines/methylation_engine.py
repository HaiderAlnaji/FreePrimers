"""
Methylation primer design engine.

Covers three bisulfite-based assay types, which have DIFFERENT and partly
OPPOSITE primer constraints:

  MSP  (Methylation-Specific PCR): two primer SETS, one matching the METHYLATED
       converted sequence and one matching the UNMETHYLATED converted sequence.
       Primers MUST cover CpG sites (so they discriminate methylation state) and
       SHOULD have a CpG at/near the 3' end for allele-style discrimination.
  BSP  (Bisulfite Sequencing PCR): one primer pair that amplifies the converted
       template REGARDLESS of methylation state, so primers MUST AVOID CpG sites
       (otherwise they'd bias toward one state). They should still include
       several converted non-CpG Cs (now Ts) to ensure they only bind fully
       converted DNA.
  MS-HRM (Methylation-Sensitive HRM): like BSP (CpG-free primers amplifying both
       states) but the amplicon spans CpGs so methylated vs unmethylated
       amplicons differ in GC content -> different melt temperature. Primers are
       BSP-style; the design goal additionally wants a short amplicon with
       enough CpGs inside to create a measurable Tm shift.

Key chemistry (verified):
  - Bisulfite: unmethylated C -> T (via U). Methylated C (CpG context) stays C.
  - After conversion the two strands are NO LONGER complementary; each strand is
    designed independently. We design on the top (converted sense) strand here
    and report; the bottom strand can be designed by reverse-complementing input.
  - "All C must be considered converted unless it is a CpG that may be
    methylated" is the core rule.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def revcomp(s: str) -> str:
    return "".join(COMP.get(b, "N") for b in reversed(s.upper()))


def cpg_positions(seq: str) -> list:
    seq = seq.upper()
    return [i for i in range(len(seq) - 1) if seq[i:i + 2] == "CG"]


def bisulfite_convert(seq: str, methylated_cpg: bool = False) -> str:
    """Convert genomic sequence as bisulfite treatment would.

    methylated_cpg=False: every C -> T (fully unmethylated).
    methylated_cpg=True:  C in CpG context stays C; all other C -> T.
    """
    seq = seq.upper()
    out = []
    n = len(seq)
    for i, b in enumerate(seq):
        if b == "C":
            is_cpg = (i + 1 < n and seq[i + 1] == "G")
            out.append("C" if (is_cpg and methylated_cpg) else "T")
        else:
            out.append(b)
    return "".join(out)


# --- simple Tm (nearest-neighbour-free Wallace/GC for bisulfite low-GC primers) ---
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
    if not seq:
        return 0.0
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


@dataclass
class Primer:
    seq: str
    start: int          # 0-based start on the converted template it was designed against
    end: int            # exclusive
    strand: str         # 'fwd' | 'rev'
    tm: float
    gc: float
    n_cpg: int          # CpG sites covered (in converted-template coordinates)
    role: str = ""
    note: str = ""


@dataclass
class AssayDesign:
    assay: str
    primers: list = field(default_factory=list)
    amplicon: Optional[dict] = None
    warnings: list = field(default_factory=list)
    notes: list = field(default_factory=list)


# ----------------------------------------------------------------- helpers
def _converted_cpg_indices(genomic: str) -> list:
    """CpG indices on the genomic sequence (these are where methylation lives)."""
    return cpg_positions(genomic)


def _score_primer_region(sub: str, tm_target: float) -> float:
    return abs(tm_basic(sub) - tm_target) + abs(gc_pct(sub) - 50) * 0.05


# ----------------------------------------------------------------- MSP
def design_msp(genomic: str, tm_target: float = 58.0,
               len_lo: int = 22, len_hi: int = 30) -> dict:
    """Design methylated and unmethylated primer SETS.

    MSP primers must (a) sit on converted template, (b) cover >=1 CpG, ideally
    with a CpG near the 3' end so the methylated vs unmethylated base difference
    (C vs T at the CpG-C) drives allele-style discrimination.
    Returns {'methylated': AssayDesign, 'unmethylated': AssayDesign}.
    """
    genomic = genomic.upper()
    cpgs = _converted_cpg_indices(genomic)
    out = {}
    for state, meth in (("methylated", True), ("unmethylated", False)):
        conv = bisulfite_convert(genomic, methylated_cpg=meth)
        design = AssayDesign(assay=f"MSP-{state}")
        # forward primer: pick a window containing a CpG with the CpG-C as close
        # to the 3' end as a CpG allows.
        fwd = _best_msp_primer(conv, genomic, cpgs, "fwd", tm_target, len_lo, len_hi)
        rev = _best_msp_primer(conv, genomic, cpgs, "rev", tm_target, len_lo, len_hi)
        if fwd:
            fwd.role = f"{state} forward"
            design.primers.append(fwd)
        if rev:
            rev.role = f"{state} reverse"
            design.primers.append(rev)
        if not cpgs:
            design.warnings.append("No CpG sites in the region — MSP is not "
                                   "applicable (nothing to discriminate).")
        if fwd and rev:
            amp_lo, amp_hi = fwd.start, rev.end
            design.amplicon = {"start": amp_lo, "end": amp_hi, "size": amp_hi - amp_lo}
        design.notes.append("MSP primers cover CpG sites; the methylated set "
                            "matches C at CpGs, the unmethylated set matches T. "
                            "Run both sets on bisulfite-converted DNA; the set "
                            "that amplifies reports the methylation state.")
        out[state] = design
    return out


def _best_msp_primer(conv, genomic, cpgs, strand, tm_target, len_lo, len_hi):
    """Pick an MSP primer covering a CpG, CpG-C as near the 3' end as possible."""
    n = len(conv)
    best, best_sc = None, 1e9
    if strand == "fwd":
        # forward: 3' end should be at or just after a CpG-C in the 5' half
        for cpg in cpgs:
            if cpg > n * 0.6:
                continue
            for L in range(len_lo, len_hi + 1):
                end = cpg + 2  # include the CpG (C and G) near 3' end
                start = end - L
                if start < 0 or end > n:
                    continue
                sub = conv[start:end]
                covered = sum(1 for c in cpgs if start <= c < end)
                if covered == 0:
                    continue
                sc = _score_primer_region(sub, tm_target) - covered * 1.5
                if sc < best_sc:
                    best_sc = sc
                    best = Primer(sub, start, end, "fwd", round(tm_basic(sub), 1),
                                  round(gc_pct(sub), 0), covered)
    else:
        # reverse: design on revcomp of converted template, 3' near a 3'-side CpG
        for cpg in cpgs:
            if cpg < n * 0.4:
                continue
            for L in range(len_lo, len_hi + 1):
                start = cpg
                end = start + L
                if end > n:
                    continue
                sub = conv[start:end]
                covered = sum(1 for c in cpgs if start <= c < end)
                if covered == 0:
                    continue
                rc = revcomp(sub)
                sc = _score_primer_region(rc, tm_target) - covered * 1.5
                if sc < best_sc:
                    best_sc = sc
                    best = Primer(rc, start, end, "rev", round(tm_basic(rc), 1),
                                  round(gc_pct(rc), 0), covered)
    return best


# ----------------------------------------------------------------- BSP
def design_bsp(genomic: str, tm_target: float = 55.0,
               len_lo: int = 20, len_hi: int = 30) -> AssayDesign:
    """One primer pair amplifying converted DNA regardless of methylation.

    Constraint OPPOSITE to MSP: primers must AVOID CpG sites (so they bind both
    methylated and unmethylated equally). Design against the unmethylated-
    converted template (worst case for Tm, all C->T) and require zero CpGs under
    each primer.
    """
    genomic = genomic.upper()
    conv = bisulfite_convert(genomic, methylated_cpg=False)  # all C->T
    cpgs = set(_converted_cpg_indices(genomic))
    design = AssayDesign(assay="BSP")
    n = len(conv)

    def cpg_free(start, end):
        return not any(start <= c < end for c in cpgs)

    fwd = _best_cpgfree_primer(conv, cpgs, "fwd", tm_target, len_lo, len_hi, region=(0, int(n * 0.5)))
    rev = _best_cpgfree_primer(conv, cpgs, "rev", tm_target, len_lo, len_hi, region=(int(n * 0.5), n))
    if fwd:
        fwd.role = "BSP forward"; design.primers.append(fwd)
    if rev:
        rev.role = "BSP reverse"; design.primers.append(rev)
    if not fwd or not rev:
        # diagnose: is the region too CpG-dense for CpG-free primers?
        longest = 0; cur = 0
        for i in range(n):
            if i in cpgs or (i - 1) in cpgs:
                cur = 0
            else:
                cur += 1; longest = max(longest, cur)
        design.warnings.append(
            f"Could not place CpG-free primer(s) ({'forward ' if not fwd else ''}"
            f"{'reverse' if not rev else ''}). Longest CpG-free run is {longest} nt "
            f"(need >={len_lo}). This region is too CpG-dense for classic BSP — "
            "extend the region to include CpG-free flanks, or use MSP/MS-HRM instead.")
    if fwd and rev:
        size = rev.end - fwd.start
        cpgs_in_amp = sum(1 for c in cpgs if fwd.start <= c < rev.end)
        design.amplicon = {"start": fwd.start, "end": rev.end, "size": size,
                           "cpgs_in_amplicon": cpgs_in_amp}
        design.notes.append(f"Amplicon spans {cpgs_in_amp} CpG(s) for downstream "
                            "sequencing readout.")
        if cpgs_in_amp == 0:
            design.warnings.append("Amplicon contains no CpGs — nothing to "
                                   "sequence for methylation. Choose a region "
                                   "with internal CpGs.")
    design.notes.append("BSP primers avoid CpG sites so they amplify methylated "
                        "and unmethylated templates equally; methylation is read "
                        "out by sequencing the amplicon.")
    return design


def _best_cpgfree_primer(conv, cpgs, strand, tm_target, len_lo, len_hi, region):
    n = len(conv)
    r_lo, r_hi = region
    best, best_sc = None, 1e9
    for start in range(r_lo, min(r_hi, n - len_lo)):
        for L in range(len_lo, len_hi + 1):
            end = start + L
            if end > n:
                continue
            if any(start <= c < end for c in cpgs):
                continue  # must be CpG-free
            sub = conv[start:end]
            # require the primer to contain at least a couple of converted Cs
            # (now Ts) downstream of a former C, ensuring conversion specificity:
            seqv = sub if strand == "fwd" else revcomp(sub)
            sc = _score_primer_region(seqv, tm_target)
            # prefer a non-T 3' end (T-rich 3' from conversion is less specific)
            if seqv[-1] == "T":
                sc += 0.5
            if sc < best_sc:
                best_sc = sc
                best = Primer(seqv, start, end, strand, round(tm_basic(seqv), 1),
                              round(gc_pct(seqv), 0), 0)
    return best


# ----------------------------------------------------------------- MS-HRM
def design_mshrm(genomic: str, tm_target: float = 55.0,
                 len_lo: int = 20, len_hi: int = 28,
                 amp_lo: int = 60, amp_hi: int = 150) -> AssayDesign:
    """MS-HRM: BSP-style CpG-free primers, but optimised for a SHORT amplicon
    that contains several CpGs, so methylated vs unmethylated amplicons differ
    enough in GC/Tm for HRM to resolve.
    """
    genomic = genomic.upper()
    conv = bisulfite_convert(genomic, methylated_cpg=False)
    cpgs = set(_converted_cpg_indices(genomic))
    n = len(conv)
    design = AssayDesign(assay="MS-HRM")

    # find a CpG-rich short window, then place CpG-free primers flanking it
    best = None
    best_sc = -1
    for s in range(0, n - amp_lo):
        for size in range(amp_lo, amp_hi + 1):
            e = s + size
            if e > n:
                break
            cpgs_in = sum(1 for c in cpgs if s <= c < e)
            if cpgs_in == 0:
                continue
            density = cpgs_in / size
            if density > best_sc:
                best_sc = density
                best = (s, e, cpgs_in)
    if not best:
        design.warnings.append("No CpG-containing window found for an HRM amplicon.")
        return design
    s, e, cpgs_in = best
    fwd = _best_cpgfree_primer(conv, cpgs, "fwd", tm_target, len_lo, len_hi,
                               region=(max(0, s - 30), s))
    rev = _best_cpgfree_primer(conv, cpgs, "rev", tm_target, len_lo, len_hi,
                               region=(e, min(n, e + 30)))
    if fwd:
        fwd.role = "MS-HRM forward"; design.primers.append(fwd)
    if rev:
        rev.role = "MS-HRM reverse"; design.primers.append(rev)
    if not fwd or not rev:
        design.warnings.append(
            "Could not place CpG-free flanking primers around the CpG-rich window. "
            "MS-HRM needs CpG-free primer sites flanking the CpG island — provide a "
            "region that includes CpG-free flanking sequence on both sides of the island.")
    if fwd and rev:
        size = rev.end - fwd.start
        # estimate Tm difference between fully methylated and unmethylated amplicons
        amp_unmeth = bisulfite_convert(genomic[fwd.start:rev.end], methylated_cpg=False)
        amp_meth = bisulfite_convert(genomic[fwd.start:rev.end], methylated_cpg=True)
        dtm = tm_basic(amp_meth) - tm_basic(amp_unmeth)
        design.amplicon = {"start": fwd.start, "end": rev.end, "size": size,
                           "cpgs_in_amplicon": sum(1 for c in cpgs if fwd.start <= c < rev.end),
                           "est_tm_shift_meth_vs_unmeth": round(dtm, 1)}
        design.notes.append(f"Estimated melt-temperature shift between fully "
                            f"methylated and unmethylated amplicon: ~{round(dtm,1)}C "
                            "(methylated amplicon retains more C/G, melts higher).")
        if abs(dtm) < 1.0:
            design.warnings.append("Predicted Tm shift <1C — HRM may not resolve "
                                   "the two states well; pick a more CpG-dense region.")
    design.notes.append("MS-HRM uses CpG-free primers (amplify both states); the "
                        "amplicon's CpGs make methylated DNA melt at a higher "
                        "temperature, resolved by HRM.")
    return design


# ----------------------------------------------------------------- CpG island detection
def cpg_island_track(seq: str, window: int = 200, step: int = 1):
    """Sliding-window CpG-island scoring (Gardiner-Garden & Frommer 1987 criteria):
    a window qualifies as CpG-island-like if GC% > 50 AND observed/expected CpG
    ratio > 0.6. obs/exp = (#CpG * N) / (#C * #G) for the window.

    Returns a per-position track (gc%, oe ratio, island bool) and merged island
    intervals, for plotting a predicted CpG-island figure.
    """
    seq = seq.upper()
    n = len(seq)
    if n < window:
        window = n
    gc_track, oe_track, island_flag = [], [], []
    # precompute prefix sums for speed
    pc = [0] * (n + 1)
    pg = [0] * (n + 1)
    pcg = [0] * (n + 1)
    for i in range(n):
        pc[i + 1] = pc[i] + (1 if seq[i] == "C" else 0)
        pg[i + 1] = pg[i] + (1 if seq[i] == "G" else 0)
        pcg[i + 1] = pcg[i] + (1 if (i + 1 < n and seq[i:i + 2] == "CG") else 0)
    centers = []
    for s in range(0, n - window + 1, step):
        e = s + window
        cC = pc[e] - pc[s]
        cG = pg[e] - pg[s]
        cCG = pcg[e] - pcg[s]
        gc = 100.0 * (cC + cG) / window
        exp = (cC * cG) / window if (cC and cG) else 0.0
        oe = (cCG / exp) if exp > 0 else 0.0
        is_island = (gc > 50.0 and oe > 0.6)
        centers.append(s + window // 2)
        gc_track.append(round(gc, 1))
        oe_track.append(round(oe, 2))
        island_flag.append(is_island)
    # merge contiguous island windows into intervals (in window-center coords -> expand)
    intervals = []
    run_start = None
    for i, flag in enumerate(island_flag):
        if flag and run_start is None:
            run_start = i
        elif not flag and run_start is not None:
            intervals.append((max(0, centers[run_start] - window // 2),
                              min(n, centers[i - 1] + window // 2)))
            run_start = None
    if run_start is not None:
        intervals.append((max(0, centers[run_start] - window // 2),
                          min(n, centers[-1] + window // 2)))
    return {"centers": centers, "gc": gc_track, "oe": oe_track,
            "island": island_flag, "intervals": intervals,
            "window": window, "region_length": n}
