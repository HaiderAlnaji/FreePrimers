"""
Empirical allele-discrimination model for tetra-primer ARMS.

Grounded in measured Taq extension efficiencies, not bulk duplex dG
(which we showed gives the scientifically wrong answer). Sources:

  Huang MM, Arnheim N, Goodman MF. Extension of base mispairs by Taq
  DNA polymerase: implications for single nucleotide discrimination in
  PCR. Nucleic Acids Res. 1992;20(17):4567-4573.
    -> relative extension efficiency for all 12 3'-terminal mispairs.

  Simsek M (Simsek & Adnan). Discrimination of primer 3'-nucleotide
  mismatch by Taq DNA polymerase during PCR. (2000).
    -> 3'-terminal T/G/C give 40-100x discrimination, A is poor;
       penultimate (-2) position ~1/5 as much (8-20x); efficiency is
       reduced when T or especially A occupy the penultimate position.

  Multiple sources (e.g. Lefever; Stadhouders; recent polymerase
  mismatch studies): >1 contiguous 3'-terminal mismatch causes near-
  complete amplification blockage. This is the mechanistic basis for
  the deliberate -2 mismatch in ARMS.

[VERIFY] All numeric values are transcribed from the above for
benchmarking; confirm against the primary sources before any appear in
a manuscript.

NOTE: values are *orders of magnitude* (log10 relative extension
efficiency). They encode the documented ranking, not lab-exact rates,
and are intended to drive design *ranking* and *warnings*, then to be
checked against the published fold-discrimination data.
"""

from __future__ import annotations

from dataclasses import dataclass, field

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def comp(b: str) -> str:
    return COMP[b]


# --- Huang 1992: -log10(relative extension efficiency) for primer.template
#     mispairs. Higher = more strongly suppressed = better discrimination.
#     A correct pair = 0. Keys are (primer_base, template_base).
HUANG_NEG_LOG10 = {
    # transitions: 10^-3 .. 10^-4  -> ~3.5
    ("A", "C"): 3.5, ("C", "A"): 3.5, ("G", "T"): 3.5, ("T", "G"): 3.5,
    # transversions
    ("T", "C"): 4.5, ("T", "T"): 4.5,          # 10^-4 .. 10^-5
    ("A", "A"): 6.0,                              # ~10^-6
    ("A", "G"): 6.5, ("G", "A"): 6.5, ("G", "G"): 6.5, ("C", "C"): 6.5,  # <10^-6
    ("C", "T"): 2.0,                              # ~10^-2  (the easy exception!)
}


def mispair_suppression(primer_base: str, template_base: str) -> float:
    """-log10 relative extension efficiency for a primer.template 3' pair.

    0.0 for a Watson-Crick match; larger for more strongly suppressed
    mispairs. Returns 0.0 if the pair is complementary (a match).
    """
    if template_base == comp(primer_base):
        return 0.0
    return HUANG_NEG_LOG10[(primer_base, template_base)]


# --- Simsek & Adnan: penultimate (-2) deliberate-mismatch contribution.
#     ~1/5 of terminal discrimination (8-20x ~ 0.9-1.3 log). We use ~1.0
#     log of *added* suppression from a -2 mismatch on the wrong allele,
#     plus a "double-mismatch blockage" boost when the wrong allele ends
#     up with both the terminal AND the -2 mismatched (near-total block).
PENULT_ADDED_LOG = 1.0
# A contiguous double 3' mismatch on the wrong allele adds extra blockage, but
# kept SMALL so the well-measured terminal-mispair value (Huang) remains the
# dominant, differentiating signal rather than saturating every design.
DOUBLE_MISMATCH_BLOCK_BONUS = 0.5

# Matched-allele extension penalties (Simsek & Adnan): a -2 deliberate
# mismatch reduces the on-target efficiency too; T and especially A at
# the -2 position are worse. Expressed as log10 efficiency LOST on the
# matched template (subtracted from a matched-primer "budget").
PENULT_MATCHED_PENALTY = {"A": 1.3, "T": 1.0, "C": 0.6, "G": 0.6}
# A matched primer should retain enough extension to prime; require the
# matched-allele net log-efficiency loss to stay within this budget.
MATCHED_LOSS_BUDGET = 1.6


@dataclass
class AlleleDesign:
    detected: str                 # sense allele this primer reports
    other: str                    # the competing sense allele
    orientation: str              # 'forward' (primer == sense) or 'reverse'
    deliberate_minus2: str        # base placed at -2 (deliberate mismatch)
    discrimination_log: float = 0.0   # net log10 fold suppression of wrong allele
    matched_loss_log: float = 0.0     # log10 efficiency lost on the matched allele
    primable: bool = True
    warnings: list = field(default_factory=list)
    detail: dict = field(default_factory=dict)


def _terminal_primer_and_template(detected: str, other: str, orientation: str):
    """Return (primer_terminal_base, template_base_matched, template_base_wrong).

    forward primer == sense: 3' base = detected; anneals to anti-sense,
      so template base = comp(allele_present).
    reverse primer == anti-sense: 3' base = comp(detected); anneals to
      sense, so template base = allele_present (sense).
    """
    if orientation == "forward":
        p_term = detected
        t_matched = comp(detected)
        t_wrong = comp(other)
    else:  # reverse
        p_term = comp(detected)
        t_matched = detected
        t_wrong = other
    return p_term, t_matched, t_wrong


def evaluate(detected: str, other: str, orientation: str, deliberate_minus2: str,
             genomic_minus2_sense: str, rt_rna_template: bool = False) -> AlleleDesign:
    """Score one allele-specific primer design.

    genomic_minus2_sense: the true sense-strand base immediately 5' of the
      SNP (the -2 position on a forward primer). The deliberate base must
      differ from the base that would correctly pair there, or it is not a
      real mismatch.
    rt_rna_template: if True, the discriminating mismatch forms against an RNA
      template during reverse transcription (miRNA context). A primer-G vs
      RNA-U is then a G.U WOBBLE pair, which RTs extend far more readily than
      the corresponding DNA G.T mispair (Chen et al. 2005: let-7 G-U
      cross-reaction during RT). We down-weight terminal G.(U) suppression in
      this context. For DNA qPCR templates (default) the Huang value stands.
    """
    d = AlleleDesign(detected=detected, other=other, orientation=orientation,
                     deliberate_minus2=deliberate_minus2)

    p_term, t_matched, t_wrong = _terminal_primer_and_template(detected, other, orientation)

    # --- terminal mismatch suppression on the WRONG allele (Huang) ---
    term_sup = mispair_suppression(p_term, t_wrong)

    # RNA-template (RT-step) wobble correction: primer G : template U(=T) is a
    # G.U wobble that barely destabilises and is poorly discriminated. Cap it.
    if rt_rna_template and (p_term, t_wrong) == ("G", "T"):
        term_sup = min(term_sup, 1.0)
        d.warnings.append("Discriminating mismatch is a 3'-terminal G.U wobble on an "
                          "RNA template (RT step); wobble pairs extend readily and "
                          "discriminate poorly (Chen et al. 2005). Expect residual "
                          "cross-reactivity; consider an orthogonal approach (LNA).")
    d.detail["terminal_mispair"] = f"{p_term}.{t_wrong}"
    d.detail["terminal_suppression_log"] = term_sup

    # Terminal A is a poor discriminator (Simsek & Adnan): cap its benefit.
    if p_term == "A":
        d.warnings.append("3'-terminal A is a weak discriminator (Simsek & Adnan); "
                          "consider the opposite strand for this allele.")
    # C.T exception (Huang): primer C vs template T leaks ~10^-2.
    if (p_term, t_wrong) == ("C", "T"):
        d.warnings.append("Wrong-allele mispair is C.T, which Taq extends readily "
                          "(~10^-2); poor discrimination on this strand.")

    # --- deliberate -2 mismatch contribution ---
    # On the matched allele the -2 base must mismatch the template to be a
    # 'deliberate mismatch'. The template -2 base (for a forward primer) is
    # comp(genomic_minus2_sense).
    if orientation == "forward":
        template_minus2 = comp(genomic_minus2_sense)
    else:
        template_minus2 = genomic_minus2_sense  # reverse primer reads sense strand
    deliberate_is_mismatch = (template_minus2 != comp(deliberate_minus2))

    added = 0.0
    if deliberate_is_mismatch:
        added += PENULT_ADDED_LOG
        # wrong allele now has BOTH terminal (mismatch) and -2 (mismatch) ->
        # contiguous double 3' mismatch -> near-complete blockage.
        if term_sup > 0:
            added += DOUBLE_MISMATCH_BLOCK_BONUS
        d.matched_loss_log = PENULT_MATCHED_PENALTY.get(deliberate_minus2, 0.8)
    else:
        d.warnings.append("Chosen -2 base matches the template; it is NOT a "
                          "deliberate mismatch and adds no ARMS discrimination.")

    d.discrimination_log = round(term_sup + added, 2)
    d.matched_loss_log = round(d.matched_loss_log, 2)
    d.primable = d.matched_loss_log <= MATCHED_LOSS_BUDGET
    if not d.primable:
        d.warnings.append(f"Matched-allele efficiency loss ({d.matched_loss_log} log) "
                          f"exceeds budget ({MATCHED_LOSS_BUDGET}); matched primer may "
                          f"prime poorly. Prefer a milder -2 base (C/G).")
    return d


def best_design(detected: str, other: str, genomic_minus2_sense: str):
    """Search orientation x deliberate-(-2) base for the best discrimination.

    Returns (best AlleleDesign, [all primable AlleleDesigns sorted]).
    Objective: maximise wrong-allele suppression subject to the matched
    allele remaining primable.
    """
    cands = []
    for orientation in ("forward", "reverse"):
        if orientation == "forward":
            template_minus2 = comp(genomic_minus2_sense)
        else:
            template_minus2 = genomic_minus2_sense
        for d2 in "ACGT":
            # require it to be a real mismatch at -2
            if template_minus2 == comp(d2):
                continue
            cands.append(evaluate(detected, other, orientation, d2, genomic_minus2_sense))
    primable = [c for c in cands if c.primable]
    pool = primable if primable else cands
    pool_sorted = sorted(pool, key=lambda c: (c.discrimination_log, -c.matched_loss_log), reverse=True)
    return pool_sorted[0], pool_sorted
