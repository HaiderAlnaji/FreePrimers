from __future__ import annotations

from fastapi import APIRouter

from app.models.schemas import (
    DuplexRequest,
    DuplexResult,
    HairpinRequest,
    HairpinResult,
    PanelScreenRequest,
    RnaFoldRequest,
    RnaFoldResult,
    TmRequest,
    TmResult,
)
from app.services import thermo

router = APIRouter(prefix="/thermo", tags=["thermodynamics"])


@router.post("/tm", response_model=TmResult)
def tm(req: TmRequest) -> TmResult:
    """Real nearest-neighbour Tm (primer3-py / SantaLucia core)."""
    return thermo.calc_tm(req.seq, req.salts)


@router.post("/hairpin", response_model=HairpinResult)
def hairpin(req: HairpinRequest) -> HairpinResult:
    """Real hairpin free energy with explicit structure-finding."""
    return thermo.calc_hairpin(req.seq, req.salts)


@router.post("/duplex", response_model=DuplexResult)
def duplex(req: DuplexRequest) -> DuplexResult:
    """
    Real heterodimer free energy between two oligos. Pass the same
    sequence twice for a self-dimer check.
    """
    return thermo.calc_duplex(req.seq1, req.seq2, req.salts)


@router.post("/rna-fold", response_model=RnaFoldResult)
def rna_fold(req: RnaFoldRequest) -> RnaFoldResult:
    """
    ViennaRNA MFE secondary structure for an RNA sequence — e.g. to
    check whether a pre-miRNA/stem-loop region folds in a way that
    could occlude a primer or probe binding site.
    """
    return thermo.fold_rna(req.seq)


@router.post("/panel-screen", response_model=list[DuplexResult])
def panel_screen(req: PanelScreenRequest) -> list[DuplexResult]:
    """
    All-pairs cross-reactivity screen across a named set of oligos —
    use this for an ARMS-PCR panel (all four primers together) or any
    multiplex reaction, instead of only checking each primer's
    self-dimer in isolation.
    """
    return thermo.screen_all_pairs(req.seqs, req.salts)
