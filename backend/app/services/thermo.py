"""
Thermodynamics service.

Replaces the heuristic Tm / hairpin / dimer estimates used in the
client-side FreePrimers tool with real nearest-neighbour calculations:

- primer3-py: the same C core used by Primer3 (and licensed/reused by
  most professional primer tools). Used for Tm, hairpin dG, and
  heterodimer (cross-)dG with proper structure-finding, not a
  complementarity heuristic.
- ViennaRNA: used specifically for RNA secondary structure questions
  that primer3 doesn't address — e.g. whether the *RNA* template
  (pre-miRNA / stem-loop region) folds in a way that could block
  primer or probe binding. This is a DNA-oligo tool, so ViennaRNA is
  the minority partner here; primer3 is authoritative for primer
  thermodynamics.

All energies are normalised to kcal/mol (primer3 returns cal/mol).
"""

from __future__ import annotations

from typing import Optional

import primer3
import RNA

from app.models.schemas import (
    DuplexResult,
    HairpinResult,
    RnaFoldResult,
    SaltConditions,
    TmResult,
)


def _cal_to_kcal(value: float) -> float:
    return value / 1000.0


def calc_tm(seq: str, salts: SaltConditions) -> TmResult:
    """Real SantaLucia-derived Tm via primer3's NN implementation."""
    tm = primer3.calc_tm(
        seq,
        mv_conc=salts.na_mm,
        dv_conc=salts.mg_mm,
        dntp_conc=salts.dntp_mm,
        dna_conc=salts.primer_nm,
    )
    return TmResult(seq=seq, tm_c=round(tm, 2))


def calc_hairpin(seq: str, salts: SaltConditions) -> HairpinResult:
    """
    Real hairpin free energy with explicit secondary-structure
    finding (loop position, stem length), replacing the sliding-
    window complementarity heuristic used client-side.
    """
    res = primer3.calc_hairpin(
        seq,
        mv_conc=salts.na_mm,
        dv_conc=salts.mg_mm,
        dntp_conc=salts.dntp_mm,
        dna_conc=salts.primer_nm,
    )
    return HairpinResult(
        seq=seq,
        structure_found=res.structure_found,
        dg_kcal_mol=round(_cal_to_kcal(res.dg), 2) if res.structure_found else 0.0,
        tm_c=round(res.tm, 2) if res.structure_found else None,
    )


def calc_duplex(seq1: str, seq2: str, salts: SaltConditions) -> DuplexResult:
    """
    Real heterodimer (cross-hybridisation) free energy between two
    oligos. With seq1 == seq2 this is the self-dimer check; with
    seq1 != seq2 it checks cross-reactivity between any two primers
    in a reaction (e.g. the four ARMS-PCR primers against each
    other, or forward vs. reverse in a standard pair).
    """
    res = primer3.calc_heterodimer(
        seq1,
        seq2,
        mv_conc=salts.na_mm,
        dv_conc=salts.mg_mm,
        dntp_conc=salts.dntp_mm,
        dna_conc=salts.primer_nm,
    )
    return DuplexResult(
        seq1=seq1,
        seq2=seq2,
        structure_found=res.structure_found,
        dg_kcal_mol=round(_cal_to_kcal(res.dg), 2) if res.structure_found else 0.0,
        tm_c=round(res.tm, 2) if res.structure_found else None,
    )


def fold_rna(seq: str) -> RnaFoldResult:
    """
    ViennaRNA minimum-free-energy secondary structure for an RNA
    sequence (e.g. the pre-miRNA hairpin region around a primer
    binding site). Returns dot-bracket notation and MFE in kcal/mol.
    Input is treated as RNA: T is read as U.
    """
    rna_seq = seq.upper().replace("T", "U")
    fc = RNA.fold_compound(rna_seq)
    structure, mfe = fc.mfe()
    return RnaFoldResult(seq=rna_seq, structure=structure, mfe_kcal_mol=round(mfe, 2))


def screen_all_pairs(seqs: dict[str, str], salts: SaltConditions) -> list[DuplexResult]:
    """
    All-pairs duplex screen across a set of named oligos (e.g. the
    four ARMS-PCR primers, or a primer panel for multiplex PCR).
    Includes self-dimers (i == j) and every cross-pair once.
    """
    names = list(seqs.keys())
    results: list[DuplexResult] = []
    for i, name_a in enumerate(names):
        for name_b in names[i:]:
            res = calc_duplex(seqs[name_a], seqs[name_b], salts)
            res.label = f"{name_a} \u00d7 {name_b}" if name_a != name_b else f"{name_a} (self)"
            results.append(res)
    return results
