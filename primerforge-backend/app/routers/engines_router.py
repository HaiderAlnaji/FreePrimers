"""
Engine endpoints: the three primer-design engines (ARMS, miRNA stem-loop,
SNP-aware placement), plus rsID fetch, FASTA parsing, and a validation panel
that scores any pasted primer set (ours, Primer1, Primer3...) on identical
metrics for honest head-to-head comparison.
"""

from __future__ import annotations

from fastapi import APIRouter

from app.engines import (arms_optimizer as arms, stemloop_assay,
                         snp_placement as snp, empirical_discrimination as emp)
from app.services import snp_fetch
from app.models.schemas import (
    ArmsDesignRequest, ArmsDesignResult, PrimerOut,
    MirnaDesignRequest, MirnaDesignResult, MirnaSibling,
    SnpScanRequest, SnpScanResult, SnpPlacement,
    RsidFetchRequest, RsidFetchResult,
    FastaParseRequest, FastaParseResult,
    ArmsComparisonRequest, ArmsComparisonResult, ComparisonRow,
    MirnaCompareRequest, MirnaCompareResult, MirnaCompareRow,
    MirnaCompareSiblingScore,
)

router = APIRouter(prefix="/engines", tags=["engines"])

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


# ----------------------------------------------------------------- ARMS
@router.post("/arms/design", response_model=ArmsDesignResult)
def arms_design(req: ArmsDesignRequest) -> ArmsDesignResult:
    p = arms.Params(tm_target=req.tm_target)
    designs = arms.optimize(req.sense, req.snp_index, req.allele1, req.allele2, p, n_results=1)
    if not designs:
        return ArmsDesignResult(primers=[], sizes={}, discrimination={}, score=0.0,
                                warnings=["No valid design found. Provide more flanking "
                                          "sequence or adjust Tm/size constraints."])
    d = designs[0]
    primers = []
    for nm in ("FO", "FI", "RI", "RO"):
        pr = d.primers[nm]
        primers.append(PrimerOut(
            name=nm, seq=pr.seq, role=pr.role, tm=round(pr.tm, 1),
            gc=round(pr.gc, 0), hairpin=round(pr.hairpin, 2),
            detects=pr.detail.get("detects"), minus2=pr.detail.get("minus2")))
    return ArmsDesignResult(primers=primers, sizes=d.sizes,
                            discrimination=d.discrimination,
                            score=round(d.score, 2), warnings=d.warnings)


# ----------------------------------------------------------------- miRNA
@router.post("/mirna/design", response_model=MirnaDesignResult)
def mirna_design(req: MirnaDesignRequest) -> MirnaDesignResult:
    a = stemloop_assay.design_assay(req.target_name, req.family[req.target_name],
                                    req.family, forward_tail=req.forward_tail)
    sibs = [MirnaSibling(
        name=s.name, diff_pos=s.diff_pos, from_3p=s.from_3p, regime=s.regime,
        terminal_mismatch=s.terminal_mismatch,
        discrimination_log=s.discrimination_log, note=s.note) for s in a.siblings]
    return MirnaDesignResult(
        target=a.target, rt_primer=a.rt_primer, rt_extension=a.rt_extension,
        forward_primer=a.forward_primer, forward_3p_pos=a.forward_3p_pos,
        universal_reverse=a.universal_reverse, siblings=sibs,
        warnings=a.warnings, notes=a.notes)


# ----------------------------------------------------------------- SNP-aware
@router.post("/snp/scan", response_model=SnpScanResult)
def snp_scan(req: SnpScanRequest) -> SnpScanResult:
    variants = [snp.Variant(pos=v["pos"], ref=v.get("ref", "N"), alt=v.get("alt", "N"),
                            maf=float(v.get("maf", 0.0)), rsid=v.get("rsid", ""))
                for v in req.variants]
    uw = 0.5 if req.conservative else 0.15
    best = snp.scan_region(req.sense, variants, primer_len=req.primer_len,
                           strand=req.strand, region_offset=req.region_offset,
                           top=req.top, unknown_weight=uw)
    placements = [SnpPlacement(seq=p.seq, start=p.start, end=p.end, strand=p.strand,
                               total_risk=p.total_risk, verdict=p.verdict,
                               n_overlaps=len(p.overlaps)) for p in best]
    return SnpScanResult(placements=placements,
                         mode="conservative" if req.conservative else "placement-default")


# ----------------------------------------------------------------- rsID fetch
@router.post("/snp/fetch-rsid", response_model=RsidFetchResult)
def fetch_rsid(req: RsidFetchRequest) -> RsidFetchResult:
    r = snp_fetch.fetch_snp(req.rsid, flank=req.flank)
    if not r.ok:
        return RsidFetchResult(ok=False, rsid=req.rsid,
                               note=r.error or "Lookup failed. Paste or upload the sequence instead.")
    return RsidFetchResult(ok=True, rsid=r.rsid, sense=r.sense, snp_index=r.snp_index,
                           allele1=r.allele1, allele2=r.allele2,
                           assembly=r.assembly, chromosome=r.chromosome,
                           position=r.position, gene=r.gene, consequence=r.consequence,
                           maf=r.maf, ancestral=r.ancestral, note=r.note)


# ----------------------------------------------------------------- FASTA parse
@router.post("/util/parse-fasta", response_model=FastaParseResult)
def parse_fasta(req: FastaParseRequest) -> FastaParseResult:
    header, seq = "", []
    for line in req.fasta.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq:
                break  # only first record
            header = line[1:]
            continue
        seq.append(line)
    s = "".join(seq).upper().replace("U", "T")
    s = "".join(c for c in s if c in "ACGT")
    return FastaParseResult(sequence=s, header=header, length=len(s))


# ----------------------------------------------------------------- comparison panel
def _revcomp(s: str) -> str:
    return "".join(COMP.get(b, "N") for b in reversed(s.upper()))


def _score_set(sense, S, a1, a2, fo, ro, fi, ri, fi_det, ri_det, fi_m2, ri_m2):
    """Score one four-primer set on shared neutral metrics."""
    def tm(x): return arms.tm(x)
    def gc(x): return arms.gc(x)
    def hp(x): return arms.hairpin_dg(x)
    tms = [tm(fo), tm(fi), tm(ri), tm(ro)]
    tm_spread = round(max(tms) - min(tms), 2)
    # worst cross-dimer over all pairs
    seqs = [fo, fi, ri, ro]
    worst_d = 0.0
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            worst_d = min(worst_d, arms.heterodimer_dg(seqs[i], seqs[j]))
    worst_h = min(hp(fo), hp(fi), hp(ri), hp(ro))
    # discrimination via empirical model using actual detected alleles + -2
    fi_other = a2 if fi_det == a1 else a1
    ri_other = a2 if ri_det == a1 else a1
    ev_fi = emp.evaluate(fi_det, fi_other, "forward", fi_m2 or "C", sense[S - 1])
    ev_ri = emp.evaluate(ri_det, ri_other, "reverse", ri_m2 or "C", sense[S + 1])
    min_disc = round(min(ev_fi.discrimination_log, ev_ri.discrimination_log), 2)
    # amplicon sizes by locating primers (outers exact; inners by body)
    fo_pos = sense.find(fo)
    ro_pos = sense.find(_revcomp(ro))
    ro_5 = ro_pos + len(ro) - 1 if ro_pos >= 0 else -1
    band_sep = None
    in_range = False
    if fo_pos >= 0 and ro_pos >= 0:
        size_fi = (ro_5 - S) + len(fi)
        size_ri = (S - fo_pos) + len(ri)
        band_sep = abs(size_fi - size_ri)
        in_range = (80 <= size_fi <= 400) and (80 <= size_ri <= 400)
    return tm_spread, min_disc, band_sep, round(worst_d, 2), round(worst_h, 2), in_range


@router.post("/arms/compare", response_model=ArmsComparisonResult)
def arms_compare(req: ArmsComparisonRequest) -> ArmsComparisonResult:
    rows = []
    for s in req.sets:
        try:
            ts, md, bs, wd, wh, ir = _score_set(
                req.sense, req.snp_index, req.allele1, req.allele2,
                s.fo.upper(), s.ro.upper(), s.fi.upper(), s.ri.upper(),
                s.fi_detects.upper(), s.ri_detects.upper(),
                s.fi_minus2.upper(), s.ri_minus2.upper())
            rows.append(ComparisonRow(label=s.label, tm_spread=ts, min_disc=md,
                                      band_sep=bs, worst_dimer=wd, worst_hairpin=wh,
                                      in_range=ir))
        except Exception as e:
            rows.append(ComparisonRow(label=s.label, tm_spread=0, min_disc=0,
                                      band_sep=None, worst_dimer=0, worst_hairpin=0,
                                      in_range=False, note=f"scoring error: {e}"))
    return ArmsComparisonResult(rows=rows, metric_help={
        "tm_spread": "max-min Tm across 4 primers (lower better)",
        "min_disc": "min predicted wrong-allele suppression, log-fold (higher better)",
        "band_sep": "allele band size difference in bp (higher better, >40 ideal)",
        "worst_dimer": "most negative cross-dimer dG kcal/mol (closer to 0 better)",
        "worst_hairpin": "most negative hairpin dG (closer to 0 better)",
        "in_range": "both allele amplicons within 80-400 bp",
    })


# ----------------------------------------------------------------- miRNA compare
def _to_dna(s: str) -> str:
    return s.strip().upper().replace("U", "T")


def _forward_3p_on_target(forward_primer: str, target_dna: str):
    """Find where a forward primer's 3' terminus maps onto the target miRNA.

    The forward primer may carry a 5' tail, so we match the LONGEST suffix of
    the primer that occurs in the target and ends as far 3' as possible.
    Returns the 0-based target index of the primer's 3' end, or None.
    """
    fp = _to_dna(forward_primer)
    best = None
    # try suffixes of decreasing length (>=8 to avoid spurious short matches)
    for start in range(len(fp)):
        suf = fp[start:]
        if len(suf) < 8:
            break
        idx = target_dna.find(suf)
        if idx >= 0:
            end = idx + len(suf) - 1
            if best is None or end > best:
                best = end
    return best


@router.post("/mirna/compare", response_model=MirnaCompareResult)
def mirna_compare(req: MirnaCompareRequest) -> MirnaCompareResult:
    target = _to_dna(req.family[req.target_name])
    n = len(target)
    rt_window_start = n - 6  # positions >= this are inside the RT 6-nt anneal
    sibs = {k: _to_dna(v) for k, v in req.family.items() if k != req.target_name}

    rows = []
    for cand in req.candidates:
        end3p = _forward_3p_on_target(cand.forward_primer, target)
        if end3p is None:
            rows.append(MirnaCompareRow(
                label=cand.label, note="Forward primer does not map to the target "
                "miRNA (check sequence/strand)."))
            continue
        reaches = end3p < rt_window_start
        sib_scores = []
        worst = 99.0
        for name, s in sibs.items():
            diffs = [i for i in range(min(n, len(s))) if target[i] != s[i]]
            if not diffs:
                sib_scores.append(MirnaCompareSiblingScore(
                    sibling=name, regime="identical", discrimination_log=0.0))
                worst = min(worst, 0.0)
                continue
            nearest = max(diffs)
            # discrimination depends on whether THIS primer's 3' reaches the diff
            if nearest <= end3p and nearest < rt_window_start:
                ev = emp.evaluate(target[nearest], s[nearest], "forward",
                                  target[nearest - 1] if nearest > 0 else "A",
                                  target[nearest - 1] if nearest > 0 else "A",
                                  rt_rna_template=False)
                regime = "forward-terminal"
                disc = round(ev.discrimination_log, 2)
            else:
                ev = emp.evaluate(target[nearest], s[nearest], "forward",
                                  target[nearest - 1] if nearest > 0 else "A",
                                  target[nearest - 1] if nearest > 0 else "A",
                                  rt_rna_template=True)
                regime = "rt-window"
                disc = round(min(ev.discrimination_log, 1.0), 2)
            sib_scores.append(MirnaCompareSiblingScore(
                sibling=name, regime=regime, discrimination_log=disc))
            worst = min(worst, disc)
        rows.append(MirnaCompareRow(
            label=cand.label, forward_3p_pos=end3p, reaches_terminal=reaches,
            min_discrimination=round(worst, 2) if worst < 99 else 0.0,
            siblings=sib_scores))
    return MirnaCompareResult(target=req.target_name, rows=rows)


# ----------------------------------------------------------------- miRNA fetch (bundled DB)
from app.data.mirna_db import search_mirna as _mirna_search
from app.models.schemas import MirnaFetchRequest, MirnaFetchResult, MirnaFetchRecord


@router.post("/mirna/fetch", response_model=MirnaFetchResult)
def mirna_fetch(req: MirnaFetchRequest) -> MirnaFetchResult:
    """Search the bundled miRBase 22.1 human mature miRNA database.

    Accepts miRNA name (hsa-miR-21-5p, miR-21, let-7a, let7a),
    family name (let-7, miR-155), or MIMAT accession (MIMAT0000076).
    Returns matched mature sequences ready to paste into the miRNA design panel.
    """
    hits = _mirna_search(req.query.strip(), req.max_results)
    if not hits:
        return MirnaFetchResult(
            ok=False, query=req.query, hits=[], total=0,
            note=f"No miRNA found matching '{req.query}'. Try: name (miR-21, let-7a), "
                 "family (let-7), or MIMAT accession (MIMAT0000076). "
                 "Database covers 194 human mature miRNAs (miRBase 22.1).")
    records = [MirnaFetchRecord(name=h.name, accession=h.accession, sequence=h.sequence,
                                family=h.family, arm=h.arm, note=h.note) for h in hits]
    return MirnaFetchResult(ok=True, query=req.query, hits=records, total=len(records),
                            source="miRBase 22.1 (bundled)",
                            note="Sequences sourced from miRBase 22.1 (Kozomara et al., 2019). "
                                 "Confirm against current miRBase before manuscript submission.")



from app.engines import methylation_engine as methyl, hrm_engine as hrm
from app.models.schemas import (MethylationRequest, MethylationResult, MethylPrimerOut,
                                HrmRequest, HrmResult)


def _mp(p):
    return MethylPrimerOut(seq=p.seq, start=p.start, end=p.end, strand=p.strand,
                           tm=p.tm, gc=p.gc, n_cpg=getattr(p, "n_cpg", 0),
                           role=p.role, note=p.note).model_dump()


@router.post("/methylation/design", response_model=MethylationResult)
def methylation_design(req: MethylationRequest) -> MethylationResult:
    cpgs = methyl.cpg_positions(req.sense)
    n_cpg = len(cpgs)
    conv_u = methyl.bisulfite_convert(req.sense, methylated_cpg=False)
    conv_m = methyl.bisulfite_convert(req.sense, methylated_cpg=True)
    win = 200 if len(req.sense) >= 200 else max(50, len(req.sense) // 2)
    island = methyl.cpg_island_track(req.sense, window=win, step=max(1, len(req.sense) // 300))
    extra = dict(cpg_positions=cpgs, converted_unmethylated=conv_u,
                 converted_methylated=conv_m, region_length=len(req.sense),
                 cpg_island=island)
    if req.mode == "msp":
        r = methyl.design_msp(req.sense, tm_target=req.tm_target)
        sets = {state: [_mp(p) for p in d.primers] for state, d in r.items()}
        warns = sorted({w for d in r.values() for w in d.warnings})
        notes = sorted({nt for d in r.values() for nt in d.notes})
        amp = next((d.amplicon for d in r.values() if d.amplicon), None)
        return MethylationResult(assay="MSP", sets=sets, amplicon=amp,
                                 n_cpg_region=n_cpg, warnings=warns, notes=notes, **extra)
    elif req.mode == "bsp":
        d = methyl.design_bsp(req.sense, tm_target=req.tm_target)
        return MethylationResult(assay="BSP", sets={"primary": [_mp(p) for p in d.primers]},
                                 amplicon=d.amplicon, n_cpg_region=n_cpg,
                                 warnings=d.warnings, notes=d.notes, **extra)
    else:
        d = methyl.design_mshrm(req.sense, tm_target=req.tm_target)
        # MS-HRM melt curves: methylated vs unmethylated amplicon (the readout)
        ms_curves = None
        if d.amplicon and "start" in d.amplicon:
            from app.engines import hrm_engine as _hrm
            seg = req.sense[d.amplicon["start"]:d.amplicon["end"]]
            amp_u = methyl.bisulfite_convert(seg, methylated_cpg=False)
            amp_m = methyl.bisulfite_convert(seg, methylated_cpg=True)
            tm_u, tm_m = _hrm._amplicon_tm(amp_u), _hrm._amplicon_tm(amp_m)
            tlo, thi = round(min(tm_u, tm_m) - 8), round(max(tm_u, tm_m) + 5)
            xs, mu = _hrm._melt_trace(tm_u, tlo, thi, 0.2)
            _, mm = _hrm._melt_trace(tm_m, tlo, thi, 0.2)
            ms_curves = {"temperature": xs,
                         "melt": {"unmethylated": [round(v, 4) for v in mu],
                                  "methylated": [round(v, 4) for v in mm]},
                         "difference": {"unmethylated": _hrm._deriv(xs, mu),
                                        "methylated": _hrm._deriv(xs, mm)},
                         "tm": {"unmethylated": round(tm_u, 2), "methylated": round(tm_m, 2)},
                         "model": "Predicted methylated vs unmethylated amplicon melt "
                                  "(empirical product Tm, illustrative — not instrument data). "
                                  "Methylated DNA retains the CpG C/G pairs (higher GC) and "
                                  "melts higher; the Tm gap is the MS-HRM discrimination signal."}
        return MethylationResult(assay="MS-HRM", sets={"primary": [_mp(p) for p in d.primers]},
                                 amplicon=d.amplicon, n_cpg_region=n_cpg,
                                 melt_curves=ms_curves,
                                 warnings=d.warnings, notes=d.notes, **extra)


@router.post("/hrm/design", response_model=HrmResult)
def hrm_design(req: HrmRequest) -> HrmResult:
    if req.mode == "standard":
        d = hrm.design_standard_hrm(req.sense, req.snp_index, req.allele1, req.allele2,
                                    tm_target=req.tm_target)
    else:
        d = hrm.design_as_hrm(req.sense, req.snp_index, req.allele1, req.allele2,
                              tm_target=req.tm_target)
    # predicted melt curves (illustrative)
    curves = None
    if d.amplicon and "start" in d.amplicon and "size" in d.amplicon:
        if req.mode == "standard":
            amp_seq = req.sense[d.amplicon["start"]:d.amplicon["end"]]
            snp_off = req.snp_index - d.amplicon["start"]
            if 0 <= snp_off < len(amp_seq):
                curves = hrm.predict_melt_curves(amp_seq, snp_off, req.allele1, req.allele2)
        else:
            # single-tube AS-HRM: two allele products of different Tm (5' GC-tag shift).
            # Three genotype curves: a1/a1 (low peak), a2/a2 (high peak), a1/a2 (both).
            if "product_a1_tm" in d.amplicon and "product_a2_tm" in d.amplicon:
                curves = hrm.predict_as_hrm_curves(
                    d.amplicon["product_a1_tm"], d.amplicon["product_a2_tm"],
                    req.allele1, req.allele2)
    return HrmResult(mode=d.mode, snp_class=d.snp_class,
                     primers=[_mp(p) for p in d.primers], amplicon=d.amplicon,
                     melt_curves=curves, warnings=d.warnings, notes=d.notes)


# ----------------------------------------------------------------- gene/region fetch
from app.services import region_fetch
from app.models.schemas import RegionFetchRequest, RegionFetchResult as RegionFetchSchema


@router.post("/util/fetch-region", response_model=RegionFetchSchema)
def fetch_region_endpoint(req: RegionFetchRequest) -> RegionFetchSchema:
    r = region_fetch.fetch_region(req.query, max_len=req.max_len)
    if not r.ok:
        return RegionFetchSchema(ok=False, query=req.query, note=r.error)
    # serialize exons as plain dicts for JSON response
    exons_out = None
    if r.exons:
        exons_out = [{"exon_number": e.exon_number, "start": e.start, "end": e.end, "length": e.length}
                     for e in r.exons]
    return RegionFetchSchema(ok=True, query=r.query, sequence=r.sequence, region=r.region,
                             gene=r.gene, assembly=r.assembly, note=r.note,
                             exons=exons_out, transcript_id=r.transcript_id,
                             region_start=r.region_start)


from app.services import cdna_fetch
from app.models.schemas import (
    CdnaFetchRequest, CdnaFetchResultSchema,
)


@router.post("/util/fetch-cdna", response_model=CdnaFetchResultSchema)
def fetch_cdna_endpoint(req: CdnaFetchRequest) -> CdnaFetchResultSchema:
    """Fetch the spliced cDNA + cDNA-coordinate exon boundaries for a transcript."""
    r = cdna_fetch.fetch_cdna(req.query, species=req.species, transcript_id=req.transcript_id)
    if not r.ok:
        return CdnaFetchResultSchema(ok=False, query=req.query, note=r.error)
    exons_out = [{"exon_number": e.exon_number, "start": e.start, "end": e.end, "length": e.length}
                 for e in r.exons]
    return CdnaFetchResultSchema(
        ok=True, query=r.query, gene=r.gene, assembly=r.assembly,
        transcript_id=r.transcript_id, biotype=r.biotype, sequence=r.sequence,
        exons=exons_out, note=r.note,
    )


@router.post("/util/list-transcripts", response_model=CdnaFetchResultSchema)
def list_transcripts_endpoint(req: CdnaFetchRequest) -> CdnaFetchResultSchema:
    """List all transcripts of a gene so the user can choose one."""
    r = cdna_fetch.list_transcripts(req.query, species=req.species)
    if not r.ok:
        return CdnaFetchResultSchema(ok=False, query=req.query, note=r.error)
    tx_out = [{"transcript_id": t.transcript_id, "is_canonical": t.is_canonical,
               "biotype": t.biotype, "length": t.length, "n_exons": t.n_exons}
              for t in r.transcripts]
    return CdnaFetchResultSchema(ok=True, query=r.query, gene=r.gene,
                                 assembly=r.assembly, transcripts=tx_out)
