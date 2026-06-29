"""
Fixed-register nearest-neighbor (NN) thermodynamic engine for ARMS
allele discrimination.

Unlike primer3's calc_end_stability (which SLIDES to find the most
stable alignment and therefore ignores forced terminal/penultimate
mismatches), this computes the duplex free energy in a FIXED register
with the primer 3' end anchored on the SNP. Matched dinucleotide steps
use Watson-Crick NN parameters; mismatched steps use internal- or
terminal-mismatch NN parameters.

Parameters are taken from Biopython's MeltingTemp tables (no values are
transcribed by hand), which implement:
  DNA_NN4  : SantaLucia & Hicks 2004 (Watson-Crick).
  DNA_IMM1 : internal single mismatches and inosine
             (Allawi & SantaLucia 1997, 1998a-d; Peyret et al. 1999;
              Bommarito et al. 2000).
  DNA_TMM1 : terminal mismatches (Bommarito et al. 2000).
[VERIFY] confirm the table provenance/citations before manuscript use.

The extension-critical quantity for a polymerase is the stability of
the 3'-TERMINAL region, not whole-duplex Tm (a single terminal mismatch
barely moves whole-duplex Tm because the upstream stays paired). We
therefore report the cumulative NN dG of the terminal `window` steps.
"""

from __future__ import annotations

from Bio.SeqUtils import MeltingTemp as mt

NN = mt.DNA_NN4
IMM = mt.DNA_IMM1
TMM = mt.DNA_TMM1

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
R = 1.987
T37 = 310.15


def comp(b: str) -> str:
    return COMP[b]


def _is_wc(top2: str, bot2: str) -> bool:
    """Both positions Watson-Crick paired? top2 5'->3', bot2 3'->5'."""
    return comp(top2[0]) == bot2[0] and comp(top2[1]) == bot2[1]


def _lookup(top2: str, bot2: str, table) -> tuple | None:
    """Look up a NN step, trying the step and its strand-swapped equivalent.

    Biopython keys are 'XY/WZ' = top(5'->3') / bottom(3'->5'). NN params are
    invariant under reading the duplex from the other end, i.e.
    reverse(top)/reverse(bot). Try both orientations.
    """
    k1 = f"{top2}/{bot2}"
    if k1 in table:
        return table[k1]
    k2 = f"{top2[::-1]}/{bot2[::-1]}"
    if k2 in table:
        return table[k2]
    return None


def terminal_dG(primer_3: str, template_3to5: str, window: int = 5) -> float:
    """Cumulative NN dG (kcal/mol, 37C) of the terminal `window` steps.

    primer_3       : primer 3'-terminal region, 5'->3'. The LAST base is the
                     3' terminus (sits on the SNP).
    template_3to5  : template region aligned UNDER the primer, written 3'->5'
                     so template_3to5[i] pairs with primer_3[i].
    window         : number of 3'-terminal dinucleotide steps to sum.

    More negative dG = more stable 3' end = more extendable.
    """
    n = len(primer_3)
    assert len(template_3to5) == n, "primer and template windows must align 1:1"
    steps = []
    for i in range(n - 1):
        top2 = primer_3[i:i + 2]
        bot2 = template_3to5[i:i + 2]
        is_terminal = (i == n - 2)  # the 3'-terminal-most step
        wc = _is_wc(top2, bot2)
        if wc:
            params = _lookup(top2, bot2, NN)
        elif is_terminal:
            params = _lookup(top2, bot2, TMM) or _lookup(top2, bot2, IMM)
        else:
            params = _lookup(top2, bot2, IMM)
        if params is None:
            # Fallback: treat an unparameterised mismatched step as strongly
            # destabilising (no stacking contribution). Flagged via +0 dG.
            params = (0.0, 0.0)
        steps.append(params)
    # take the terminal `window` steps (closest to the 3' end)
    use = steps[-window:] if window < len(steps) else steps
    dH = sum(p[0] for p in use)
    dS = sum(p[1] for p in use)
    return dH - T37 * (dS / 1000.0)


def discrimination(primer_3: str, genomic_3to5_matched: str,
                   snp_pos_from_3prime: int = 0, other_allele_sense: str | None = None,
                   window: int = 5) -> dict:
    """ddG-based discrimination for an allele-specific primer.

    primer_3              : primer 3' region 5'->3' (terminus = detected allele
                            for a forward primer).
    genomic_3to5_matched  : matched-allele template under the primer, 3'->5'.
    snp_pos_from_3prime   : index of the SNP from the 3' terminus (0 = terminus).
    other_allele_sense    : the competing sense allele; used to build the
                            wrong-allele template by swapping the SNP base.

    Returns dG_matched, dG_wrong (terminal-window), and ddG = dG_wrong - dG_matched
    (more positive = wrong allele less stable = better discrimination).
    """
    dg_m = terminal_dG(primer_3, genomic_3to5_matched, window)
    # build wrong-allele template: at the SNP, template base = comp(allele)
    idx = (len(genomic_3to5_matched) - 1) - snp_pos_from_3prime
    tmpl_wrong = list(genomic_3to5_matched)
    tmpl_wrong[idx] = comp(other_allele_sense)
    dg_w = terminal_dG(primer_3, "".join(tmpl_wrong), window)
    return {
        "dG_matched": round(dg_m, 2),
        "dG_wrong": round(dg_w, 2),
        "ddG": round(dg_w - dg_m, 2),
    }
