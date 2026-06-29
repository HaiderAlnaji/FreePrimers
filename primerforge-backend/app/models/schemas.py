"""Pydantic models shared across routers and services."""

from __future__ import annotations

from typing import Literal, Optional

from pydantic import BaseModel, Field, field_validator

VALID_BASES = set("ACGTUacgtu")


def _clean_seq(v: str) -> str:
    v = v.strip().upper().replace("U", "T")
    if not v:
        raise ValueError("sequence cannot be empty")
    bad = set(v) - set("ACGT")
    if bad:
        raise ValueError(f"sequence contains non-ACGTU characters: {sorted(bad)}")
    return v


class SaltConditions(BaseModel):
    na_mm: float = Field(50.0, ge=0, le=1000, description="Monovalent cation conc. (mM)")
    mg_mm: float = Field(3.0, ge=0, le=50, description="Mg2+ conc. (mM)")
    dntp_mm: float = Field(0.8, ge=0, le=20, description="dNTP conc. (mM)")
    primer_nm: float = Field(250.0, ge=1, le=10000, description="Primer conc. (nM)")


class TmRequest(BaseModel):
    seq: str
    salts: SaltConditions = SaltConditions()

    @field_validator("seq")
    @classmethod
    def _seq(cls, v: str) -> str:
        return _clean_seq(v)


class TmResult(BaseModel):
    seq: str
    tm_c: float


class HairpinRequest(BaseModel):
    seq: str
    salts: SaltConditions = SaltConditions()

    @field_validator("seq")
    @classmethod
    def _seq(cls, v: str) -> str:
        return _clean_seq(v)


class HairpinResult(BaseModel):
    seq: str
    structure_found: bool
    dg_kcal_mol: float
    tm_c: Optional[float] = None


class DuplexRequest(BaseModel):
    seq1: str
    seq2: str
    salts: SaltConditions = SaltConditions()

    @field_validator("seq1", "seq2")
    @classmethod
    def _seq(cls, v: str) -> str:
        return _clean_seq(v)


class DuplexResult(BaseModel):
    seq1: str
    seq2: str
    structure_found: bool
    dg_kcal_mol: float
    tm_c: Optional[float] = None
    label: Optional[str] = None


class RnaFoldRequest(BaseModel):
    seq: str

    @field_validator("seq")
    @classmethod
    def _seq(cls, v: str) -> str:
        return _clean_seq(v)


class RnaFoldResult(BaseModel):
    seq: str
    structure: str
    mfe_kcal_mol: float


class PanelScreenRequest(BaseModel):
    """Named oligos for an all-pairs cross-reactivity screen."""

    seqs: dict[str, str] = Field(..., min_length=2, max_length=12)
    salts: SaltConditions = SaltConditions()

    @field_validator("seqs")
    @classmethod
    def _seqs(cls, v: dict[str, str]) -> dict[str, str]:
        return {name: _clean_seq(seq) for name, seq in v.items()}


# ---------------------------------------------------------------------------
# Specificity (BLAST)
# ---------------------------------------------------------------------------

SpecificityBackend = Literal["local", "ncbi", "auto"]


class SpecificityRequest(BaseModel):
    seq: str
    database: str = Field(
        "mirbase_mature",
        description="Local DB name (see /specificity/databases) or, "
        "for backend='ncbi', an NCBI database name (e.g. 'core_nt').",
    )
    backend: SpecificityBackend = "auto"
    max_hits: int = Field(25, ge=1, le=200)
    max_mismatches_3prime: int = Field(
        2, ge=0, le=10,
        description="Hits with <= this many mismatches in the 3' seed "
        "(last 5 nt) are flagged as likely to amplify off-target.",
    )
    is_mirna: bool = Field(
        False,
        description="When True, hits to miRNA/ncRNA/small-RNA records are "
        "interpreted as expected (target or family) rather than off-target risk.",
    )

    @field_validator("seq")
    @classmethod
    def _seq(cls, v: str) -> str:
        return _clean_seq(v)


class SpecificityHit(BaseModel):
    subject_id: str
    subject_title: Optional[str] = None
    pct_identity: float
    align_length: int
    mismatches: int
    gaps: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    bit_score: float
    three_prime_mismatches: Optional[int] = None
    risk: Literal["high", "moderate", "low", "expected"] = "low"
    hit_context: Optional[str] = None   # "mirna_family" | "off_target" | None


class SpecificityResult(BaseModel):
    seq: str
    database: str
    backend_used: Literal["local", "ncbi"]
    total_hits: int
    hits: list[SpecificityHit]
    warnings: list[str] = []


# ============================================================
# Engine schemas: ARMS / miRNA / SNP-aware placement, rsID fetch
# ============================================================

class ArmsDesignRequest(BaseModel):
    sense: str
    snp_index: int = Field(..., ge=0, description="0-based index of the SNP base")
    allele1: str
    allele2: str
    tm_target: float = Field(60.0, ge=40, le=75)

    @field_validator("sense")
    @classmethod
    def _s(cls, v: str) -> str:
        return _clean_seq(v)

    @field_validator("allele1", "allele2")
    @classmethod
    def _a(cls, v: str) -> str:
        v = v.strip().upper().replace("U", "T")
        if v not in {"A", "C", "G", "T"}:
            raise ValueError("allele must be a single base A/C/G/T")
        return v


class PrimerOut(BaseModel):
    name: str
    seq: str
    role: str
    tm: float
    gc: float
    hairpin: float
    detects: Optional[str] = None
    minus2: Optional[str] = None


class ArmsDesignResult(BaseModel):
    primers: list[PrimerOut]
    sizes: dict
    discrimination: dict
    score: float
    warnings: list[str] = []


class MirnaDesignRequest(BaseModel):
    target_name: str
    family: dict  # {name: sequence}
    forward_tail: str = ""


class MirnaSibling(BaseModel):
    name: str
    diff_pos: Optional[int] = None
    from_3p: Optional[int] = None
    regime: str
    terminal_mismatch: Optional[str] = None
    discrimination_log: float
    note: str = ""


class MirnaDesignResult(BaseModel):
    target: str
    rt_primer: str
    rt_extension: str
    forward_primer: str
    forward_3p_pos: int
    universal_reverse: str
    siblings: list[MirnaSibling]
    warnings: list[str] = []
    notes: list[str] = []


class SnpScanRequest(BaseModel):
    sense: str
    variants: list[dict]  # [{pos, ref, alt, maf, rsid}]
    primer_len: int = Field(20, ge=15, le=35)
    strand: Literal["fwd", "rev"] = "fwd"
    region_offset: int = 0
    top: int = Field(10, ge=1, le=50)
    conservative: bool = False

    @field_validator("sense")
    @classmethod
    def _s(cls, v: str) -> str:
        return _clean_seq(v)


class SnpPlacement(BaseModel):
    seq: str
    start: int
    end: int
    strand: str
    total_risk: float
    verdict: str
    n_overlaps: int


class SnpScanResult(BaseModel):
    placements: list[SnpPlacement]
    mode: str


class RsidFetchRequest(BaseModel):
    rsid: str
    flank: int = Field(250, ge=80, le=600)


class RsidFetchResult(BaseModel):
    ok: bool
    rsid: str
    sense: Optional[str] = None
    snp_index: Optional[int] = None
    allele1: Optional[str] = None
    allele2: Optional[str] = None
    assembly: Optional[str] = None
    chromosome: Optional[str] = None
    position: Optional[int] = None
    gene: Optional[str] = None
    consequence: Optional[str] = None
    maf: Optional[float] = None
    ancestral: Optional[str] = None
    var_class: Optional[str] = None
    note: str = ""


class FastaParseRequest(BaseModel):
    fasta: str


class FastaParseResult(BaseModel):
    sequence: str
    header: str = ""
    length: int


# --- validation/comparison panel: score arbitrary primer sets on shared metrics ---
class ComparisonPrimerSet(BaseModel):
    label: str                 # e.g. "Primer1", "Primer3", "Our optimizer"
    fo: str
    ro: str
    fi: str
    ri: str
    fi_detects: str
    ri_detects: str
    fi_minus2: str = ""
    ri_minus2: str = ""


class ArmsComparisonRequest(BaseModel):
    sense: str
    snp_index: int
    allele1: str
    allele2: str
    sets: list[ComparisonPrimerSet]

    @field_validator("sense")
    @classmethod
    def _s(cls, v: str) -> str:
        return _clean_seq(v)


class ComparisonRow(BaseModel):
    label: str
    tm_spread: float
    min_disc: float
    band_sep: Optional[int] = None
    worst_dimer: float
    worst_hairpin: float
    in_range: bool
    note: str = ""


class ArmsComparisonResult(BaseModel):
    rows: list[ComparisonRow]
    metric_help: dict


# --- miRNA comparison panel: score competing forward primers vs a family ---
class MirnaForwardCandidate(BaseModel):
    label: str
    forward_primer: str


class MirnaCompareRequest(BaseModel):
    target_name: str
    family: dict
    candidates: list[MirnaForwardCandidate]


class MirnaCompareSiblingScore(BaseModel):
    sibling: str
    regime: str
    discrimination_log: float


class MirnaCompareRow(BaseModel):
    label: str
    forward_3p_pos: Optional[int] = None
    reaches_terminal: bool = False
    min_discrimination: float = 0.0
    siblings: list[MirnaCompareSiblingScore] = []
    note: str = ""


class MirnaCompareResult(BaseModel):
    target: str
    rows: list[MirnaCompareRow]


# ============================================================
# Methylation + HRM engine schemas
# ============================================================
class MethylationRequest(BaseModel):
    sense: str
    mode: Literal["msp", "bsp", "mshrm"]
    tm_target: float = 58.0

    @field_validator("sense")
    @classmethod
    def _s(cls, v: str) -> str:
        return _clean_seq(v)


class MethylPrimerOut(BaseModel):
    seq: str
    start: int
    end: int
    strand: str
    tm: float
    gc: float
    n_cpg: int
    role: str
    note: str = ""


class MethylationResult(BaseModel):
    assay: str
    sets: dict  # {set_name: [primers]} for MSP; {"primary": [...]} for BSP/MS-HRM
    amplicon: Optional[dict] = None
    n_cpg_region: int = 0
    cpg_positions: list[int] = []
    converted_unmethylated: str = ""
    converted_methylated: str = ""
    region_length: int = 0
    cpg_island: Optional[dict] = None
    melt_curves: Optional[dict] = None
    warnings: list[str] = []
    notes: list[str] = []


class HrmRequest(BaseModel):
    sense: str
    snp_index: int = Field(..., ge=0)
    allele1: str
    allele2: str
    mode: Literal["standard", "as"]
    tm_target: float = 60.0

    @field_validator("sense")
    @classmethod
    def _s(cls, v: str) -> str:
        return _clean_seq(v)

    @field_validator("allele1", "allele2")
    @classmethod
    def _a(cls, v: str) -> str:
        v = v.strip().upper().replace("U", "T")
        if v not in {"A", "C", "G", "T"}:
            raise ValueError("allele must be a single base A/C/G/T")
        return v


class HrmResult(BaseModel):
    mode: str
    snp_class: int
    primers: list[MethylPrimerOut]
    amplicon: Optional[dict] = None
    melt_curves: Optional[dict] = None
    warnings: list[str] = []
    notes: list[str] = []


class RegionFetchRequest(BaseModel):
    query: str          # gene symbol, Ensembl id, or "chr:start-end"
    max_len: int = Field(4000, ge=100, le=10000)


class ExonInfo(BaseModel):
    exon_number: int
    start: int        # 0-based offset from region start
    end: int          # exclusive
    length: int

class CdnaExonSchema(BaseModel):
    exon_number: int
    start: int        # 0-based offset in cDNA
    end: int
    length: int

class TranscriptOptionSchema(BaseModel):
    transcript_id: str
    is_canonical: bool
    biotype: str
    length: int
    n_exons: int

class CdnaFetchRequest(BaseModel):
    query: str
    species: str = "human"
    transcript_id: str = ""

class CdnaFetchResultSchema(BaseModel):
    ok: bool
    query: str = ""
    gene: str = ""
    assembly: str = ""
    transcript_id: str = ""
    biotype: str = ""
    sequence: str = ""
    exons: Optional[list] = None
    transcripts: Optional[list] = None
    note: str = ""

class RegionFetchResult(BaseModel):
    ok: bool
    query: str
    sequence: str = ""
    region: str = ""
    gene: str = ""
    assembly: str = ""
    note: str = ""
    exons: Optional[list] = None
    transcript_id: str = ""
    region_start: int = 0


class MirnaFetchRequest(BaseModel):
    query: str                  # miRNA name, family, or MIMAT accession
    max_results: int = Field(10, ge=1, le=50)

class MirnaFetchRecord(BaseModel):
    name: str
    accession: str
    sequence: str               # mature RNA 5'->3' (U)
    family: str
    arm: str
    note: str = ""

class MirnaFetchResult(BaseModel):
    ok: bool
    query: str
    hits: list[MirnaFetchRecord] = []
    total: int = 0
    source: str = "miRBase 22.1 (bundled)"
    note: str = ""
