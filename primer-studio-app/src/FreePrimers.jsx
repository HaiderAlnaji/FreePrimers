import { useState, useMemo, useEffect, useRef } from "react";
import {
  Dna, FlaskConical, Copy, Check, ChevronDown, Settings2,
  CircleCheck, TriangleAlert, CircleAlert, Beaker, Info, Sparkles, Download,
  ShieldCheck, ShieldAlert, Loader2, Plug, PlugZap, Upload, X,
  Save, FolderOpen, Trash2, Clock, Bookmark
} from "lucide-react";

/* ============================================================
   THERMODYNAMICS ENGINE
   Nearest-neighbor parameters: SantaLucia (1998) unified set.
   Salt: SantaLucia [Na+] entropy correction + von Ahsen Mg2+
   monovalent equivalence. These are the validated models that
   Primer3/IDT build on — not the Wallace 2(A+T)+4(G+C) rule.
   ============================================================ */

const NN = {
  AA: [-7.9, -22.2], AT: [-7.2, -20.4], AC: [-8.4, -22.4], AG: [-7.8, -21.0],
  TA: [-7.2, -21.3], TT: [-7.9, -22.2], TC: [-8.2, -22.2], TG: [-8.5, -22.7],
  CA: [-8.5, -22.7], CT: [-7.8, -21.0], CC: [-8.0, -19.9], CG: [-10.6, -27.2],
  GA: [-8.2, -22.2], GT: [-8.4, -22.4], GC: [-9.8, -24.4], GG: [-8.0, -19.9],
};
const R = 1.987;
const COMP = { A: "T", T: "A", C: "G", G: "C" };

const cleanDNA = (s) => (s || "").toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, "");
const revComp = (s) => s.split("").reverse().map((b) => COMP[b] || "").join("");
const gc = (s) => (s.length ? ((s.match(/[GC]/g) || []).length / s.length) * 100 : 0);

function tm(seq, p) {
  seq = cleanDNA(seq);
  if (seq.length < 2) return null;
  let dH = 0, dS = 0;
  for (let i = 0; i < seq.length - 1; i++) {
    const v = NN[seq.substr(i, 2)];
    if (!v) return null;
    dH += v[0]; dS += v[1];
  }
  const init = (b) => (b === "A" || b === "T" ? [2.3, 4.1] : [0.1, -2.8]);
  const i5 = init(seq[0]), i3 = init(seq[seq.length - 1]);
  dH += i5[0] + i3[0]; dS += i5[1] + i3[1];
  const monoEq = p.na + 120 * Math.sqrt(Math.max(0, p.mg - p.dntp)); // mM
  dS += 0.368 * (seq.length - 1) * Math.log(monoEq / 1000);
  const ct = p.primer * 1e-9;
  return (dH * 1000) / (dS + R * Math.log(ct / 4)) - 273.15;
}

const dg37 = (pair) => {
  const v = NN[pair];
  return v ? v[0] - (310.15 * v[1]) / 1000 : 0;
};

// Heuristic dimer free energy (most stable antiparallel register).
function dimerDG(s1, s2) {
  const a = s1;
  const b = s2.split("").reverse().join("");
  let best = 0, three = false;
  for (let off = -(b.length - 1); off < a.length; off++) {
    let k = 0;
    const m = [];
    for (let i = 0; i < a.length; i++) {
      const j = i - off;
      m.push(j >= 0 && j < b.length && COMP[a[i]] === b[j] ? 1 : 0);
    }
    while (k < m.length) {
      if (m[k]) {
        const start = k;
        while (k < m.length && m[k]) k++;
        if (k - start >= 2) {
          let g = 0;
          for (let q = start; q < k - 1; q++) g += dg37(a.substr(q, 2));
          if (g < best) { best = g; three = k - 1 >= a.length - 4; }
        }
      } else k++;
    }
  }
  return { dg: best, three };
}

// Heuristic hairpin free energy (foldback stem + loop penalty).
function hairpinDG(seq) {
  const n = seq.length;
  let best = 0;
  for (let i = 0; i < n; i++) {
    for (let len = 3; len <= 8 && i + len <= n; len++) {
      const target = revComp(seq.substr(i, len));
      const idx = seq.indexOf(target, i + len + 3);
      if (idx !== -1) {
        let g = 0;
        for (let q = i; q < i + len - 1; q++) g += dg37(seq.substr(q, 2));
        const loop = idx - (i + len);
        const total = g + 3.5 + 0.4 * Math.max(0, loop - 4);
        if (total < best) best = total;
      }
    }
  }
  return best;
}

function evalPrimer(seq, p) {
  const T = tm(seq, p), G = gc(seq);
  const sd = dimerDG(seq, seq), hp = hairpinDG(seq);
  const runs = (seq.match(/A{5,}|T{5,}|G{5,}|C{5,}/g) || []).length;
  const clamp = /[GC]$/.test(seq) && !/[GC]{5,}$/.test(seq);
  const issues = [];
  if (T != null && Math.abs(T - p.target) > 3) issues.push(`Tm ${T.toFixed(1)}°C vs ${p.target}°C target`);
  if (G < 40 || G > 60) issues.push(`GC ${G.toFixed(0)}% out of 40–60%`);
  if (sd.dg < -6) issues.push(`Self-dimer ΔG ${sd.dg.toFixed(1)} kcal/mol`);
  else if (sd.three && sd.dg < -3.5) issues.push(`3′-end self-dimer`);
  if (hp < -2) issues.push(`Hairpin ΔG ${hp.toFixed(1)} kcal/mol`);
  if (runs) issues.push(`Mononucleotide run ≥5 nt`);
  const quality = issues.length === 0 ? "good" : issues.length === 1 ? "ok" : "warn";
  return { seq, tm: T, gc: G, selfDimer: sd.dg, hairpin: hp, clamp, issues, quality };
}

/* ============================================================
   BACKEND CLIENT
   Optional. With no backend configured, every function below
   resolves to null and the tool runs exactly as it does
   client-only, using the heuristic values above. Once a backend
   URL is set, Tm / hairpin / dimer values are upgraded in place to
   real primer3-py nearest-neighbour numbers, and specificity
   (BLAST) becomes available.
   ============================================================ */

const LOCAL_STORAGE_DISABLED_NOTE = "in-memory only, not persisted across reloads";
let backendUrlMemory = "";

function getBackendUrl() { return backendUrlMemory; }
function setBackendUrl(url) { backendUrlMemory = (url || "").trim().replace(/\/+$/, ""); }

async function fetchWithTimeout(url, options, ms) {
  const ctrl = new AbortController();
  const id = setTimeout(() => ctrl.abort(), ms);
  try {
    return await fetch(url, { ...options, signal: ctrl.signal });
  } finally {
    clearTimeout(id);
  }
}

async function backendFetch(path, body) {
  const base = getBackendUrl();
  if (!base) return null;
  try {
    const res = await fetchWithTimeout(`${base}${path}`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    }, 20000);
    if (!res.ok) return null;
    return await res.json();
  } catch {
    return null;
  }
}

async function checkBackendHealth(url) {
  if (!url) return false;
  try {
    const res = await fetchWithTimeout(`${url.replace(/\/+$/, "")}/health`, {}, 8000);
    if (!res.ok) return false;
    const data = await res.json();
    return data?.status === "ok";
  } catch {
    return false;
  }
}

function saltsPayload(p) {
  return { na_mm: p.na, mg_mm: p.mg, dntp_mm: p.dntp, primer_nm: p.primer };
}

async function backendTm(seq, p) {
  const r = await backendFetch("/thermo/tm", { seq, salts: saltsPayload(p) });
  return r ? r.tm_c : null;
}

async function backendHairpin(seq, p) {
  const r = await backendFetch("/thermo/hairpin", { seq, salts: saltsPayload(p) });
  return r ? r.dg_kcal_mol : null;
}

async function backendDuplex(seq1, seq2, p) {
  const r = await backendFetch("/thermo/duplex", { seq1, seq2, salts: saltsPayload(p) });
  return r ? r.dg_kcal_mol : null;
}

async function backendPanelScreen(seqs, p) {
  const r = await backendFetch("/thermo/panel-screen", { seqs, salts: saltsPayload(p) });
  return Array.isArray(r) ? r : null;
}

async function backendVerify(seq, p) {
  // One pass: real Tm, hairpin, and self-dimer for a single oligo.
  const [tmC, hp, sd] = await Promise.all([
    backendTm(seq, p),
    backendHairpin(seq, p),
    backendDuplex(seq, seq, p),
  ]);
  if (tmC == null && hp == null && sd == null) return null;
  return { tm: tmC, hairpin: hp, selfDimer: sd };
}

async function backendSpecificity(seq, database, opts = {}) {
  const base = getBackendUrl();
  if (!base) return null;
  // If the user types "ncbi" or "auto" in the database box, route to that
  // backend (remote BLAST against NCBI) instead of requiring a local database.
  const dbRaw = (database || "mirbase_mature").trim();
  const dbLower = dbRaw.toLowerCase();
  let backend = opts.backend || "local";
  let dbName = dbRaw;
  if (dbLower === "ncbi" || dbLower === "auto") {
    backend = dbLower;            // remote NCBI BLAST
    dbName = "nt";                // default NCBI nucleotide database
  }
  try {
    const res = await fetchWithTimeout(`${base}/specificity`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        seq, database: dbName,
        backend, max_hits: opts.maxHits || 25,
        is_mirna: opts.isMirna || false,
      }),
    }, 90000);
    if (!res.ok) {
      const err = await res.json().catch(() => null);
      return { error: err?.detail || `Specificity check failed (${res.status})` };
    }
    return await res.json();
  } catch (e) {
    if (e?.name === "AbortError") {
      return { error: "Specificity check timed out. If you meant local BLAST, confirm the backend was started with PRIMERFORGE_BLAST_DB_DIR set to your blastdb folder." };
    }
    return { error: "Could not reach the backend for specificity checking." };
  }
}

/* ============================================================
   miRNA DESIGN
   Stem-loop scaffold follows the Chen (2005) / Varkonyi-Gasic
   (2007) universal design. Always validate empirically.
   ============================================================ */

const SL_SCAFFOLD = "GTCGTATCCAGTGCAGGGTCCGAGGTATTCGCACTGGATACGAC";
const SL_REVERSE = "GTGCAGGGTCCGAGGT";
const TAIL = ["G", "C", "G", "C", "G", "C", "G", "C"];

function withTail(core, p) {
  let tail = "";
  for (let t = 0; t <= TAIL.length; t++) {
    const T = tm(tail + core, p);
    if (T != null && T >= p.target) break;
    if (t < TAIL.length) tail = TAIL.slice(0, t + 1).join("");
  }
  return tail;
}

function designStemLoop(miR, p) {
  const m = cleanDNA(miR);
  if (m.length < 16) return { error: "Enter a mature miRNA of at least 16 nt." };
  const overhang = revComp(m.slice(-6));
  const core = m.slice(0, Math.max(14, m.length - 6));
  const tail = withTail(core, p);
  const forward = tail + core;
  return {
    method: "Stem-loop RT-qPCR",
    blurb: "A stem-loop RT primer with a 6-nt 3′ overhang extends the short miRNA, then a miRNA-specific forward primer pairs with a universal reverse. Best specificity for closely related family members.",
    oligos: [
      { label: "Stem-loop RT primer", seq: SL_SCAFFOLD + overhang, segs: [[SL_SCAFFOLD, "scaffold"], [overhang, "overhang"]], note: "3′ overhang anneals to the miRNA 3′ end. Self-folds into a stem-loop, so a linear Tm is not meaningful here." },
      { label: "Forward primer (miRNA-specific)", seq: forward, segs: [[tail, "tail"], [core, "mir"]], evalIt: true },
      { label: "Universal reverse primer", seq: SL_REVERSE, segs: [[SL_REVERSE, "scaffold"]], evalIt: true },
    ],
    map: {
      length: m.length, seqStr: m,
      features: [
        { key: "Forward primer (miRNA-specific)", start: 0, end: core.length, label: "forward core", color: "#0d9488", strand: "fwd" },
        { key: "Stem-loop RT primer", start: m.length - 6, end: m.length, label: "RT anneal (3′ overhang)", color: "#d14545", strand: "rev" },
      ],
      spans: [],
    },
    p,
  };
}

function designPolyA(miR, p) {
  const m = cleanDNA(miR);
  if (m.length < 16) return { error: "Enter a mature miRNA of at least 16 nt." };
  const tail = withTail(m, p);
  const forward = tail + m;
  return {
    method: "Poly(A) tailing + universal reverse",
    blurb: "The miRNA is polyadenylated and reverse-transcribed with a universal oligo-dT adapter; a miRNA-specific forward primer pairs with a universal reverse. Simpler and cheaper; lower specificity than stem-loop.",
    oligos: [
      { label: "Forward primer (miRNA-specific)", seq: forward, segs: [[tail, "tail"], [m, "mir"]], evalIt: true },
      { label: "RT primer", seq: "(universal oligo-dT adapter — defined by your kit)", segs: [], note: "Poly(A)-tail the RNA, then reverse-transcribe with the kit's universal adapter primer." },
      { label: "Universal reverse primer", seq: "(adapter anchor — match to your RT adapter)", segs: [], note: "The reverse primer is the universal anchor of your chosen poly(A) RT system, not a miRNA-specific sequence." },
    ],
    map: {
      length: m.length, seqStr: m,
      features: [{ key: "Forward primer (miRNA-specific)", start: 0, end: m.length, label: "forward primer (full miRNA + 5′ tail)", color: "#0d9488", strand: "fwd" }],
      spans: [],
    },
    p,
  };
}

/* ============================================================
   STANDARD qPCR / PCR PAIR DESIGN  (a compact Primer3-style pick)
   ============================================================ */

function designStandard(template, p, exonContext = null) {
  const seq = cleanDNA(template);
  if (seq.length < 60) return { error: "Paste a target of at least 60 nt." };
  const [aMin, aMax] = p.amplicon;
  const step = Math.max(1, Math.floor(seq.length / 400));

  // ── Exon constraint helpers ──────────────────────────────────────────────
  // exonContext: { exons, strategy, selectedExons }
  const exons = exonContext?.exons;
  const strategy = exonContext?.strategy || "any_exon";
  const selectedExons = exonContext?.selectedExons || [];

  // Build per-base exon map: exonOf[i] = exon_number or null (intron)
  const exonOf = new Array(seq.length).fill(null);
  if (exons) {
    exons.forEach(ex => {
      for (let i = Math.max(0, ex.start); i < Math.min(seq.length, ex.end); i++)
        exonOf[i] = ex.exon_number;
    });
  }

  // Is every base of [start, end) in an exon?
  const allExonic = (start, end) => {
    if (!exons) return true; // no constraint without exon data
    for (let i = start; i < end; i++) if (exonOf[i] === null) return false;
    return true;
  };

  // For "neighboring exons" strategy: which exon does a primer primarily belong to?
  const primerExonNum = (start, end) => {
    if (!exons) return null;
    const counts = {};
    for (let i = start; i < end; i++) {
      const e = exonOf[i];
      if (e != null) counts[e] = (counts[e] || 0) + 1;
    }
    const keys = Object.keys(counts);
    if (!keys.length) return null;
    return Number(keys.reduce((a, b) => counts[a] > counts[b] ? a : b));
  };

  // For "junction" strategy: does a primer span an exon-intron boundary?
  const spansJunction = (start, end) => {
    if (!exons) return false;
    const first = exonOf[start], last = exonOf[end - 1];
    if (first === null || last === null) return false; // must start AND end in exon
    // Check if there's a transition (different exon numbers = crosses junction)
    for (let i = start + 1; i < end; i++)
      if (exonOf[i] !== exonOf[i - 1]) return true;
    return false;
  };

  // Should this primer pass the exon filter?
  const primerOk = (start, end, role) => {
    if (!exons || strategy === "any_exon") return allExonic(start, end);
    if (strategy === "junction") {
      return allExonic(start, end) || spansJunction(start, end);
    }
    if (strategy === "neighboring") {
      if (!allExonic(start, end)) return false;
      if (selectedExons.length === 2) {
        const myExon = primerExonNum(start, end);
        // fwd must be on selectedExons[0], rev on selectedExons[1]
        return role === "fwd" ? myExon === selectedExons[0] : myExon === selectedExons[1];
      }
      return true;
    }
    return allExonic(start, end);
  };

  const fwd = [], rev = [];
  for (let s = 0; s < seq.length - 18; s += step) {
    for (let l = 18; l <= 25 && s + l <= seq.length; l++) {
      // Forward primer: binds left-to-right at position [s, s+l)
      if (primerOk(s, s + l, "fwd")) {
        const e = evalPrimer(seq.substr(s, l), p);
        if (e.tm != null && Math.abs(e.tm - p.target) <= 3 && e.gc >= 40 && e.gc <= 60 && e.issues.length <= 1)
          fwd.push({ ...e, start: s, end: s + l, exonNum: exons ? primerExonNum(s, s + l) : null });
      }
      // Reverse primer: binds at position [s, s+l) but reads right-to-left
      // The sequence it reads is the reverse complement of seq[s..s+l]
      if (primerOk(s, s + l, "rev")) {
        const r = evalPrimer(revComp(seq.substr(s, l)), p);
        if (r.tm != null && Math.abs(r.tm - p.target) <= 3 && r.gc >= 40 && r.gc <= 60 && r.issues.length <= 1)
          rev.push({ ...r, seq: revComp(seq.substr(s, l)), bindStart: s, bindEnd: s + l,
                    exonNum: exons ? primerExonNum(s, s + l) : null });
      }
    }
  }
  const rank = (x) => x.issues.length + Math.abs(x.tm - p.target);
  fwd.sort((a, b) => rank(a) - rank(b));
  rev.sort((a, b) => rank(a) - rank(b));
  const F = fwd.slice(0, 120), V = rev.slice(0, 120);
  const pairs = [];
  for (const f of F) for (const r of V) {
    if (r.bindStart < f.end) continue;
    const amp = r.bindEnd - f.start;
    if (amp < aMin || amp > aMax) continue;
    const tmDiff = Math.abs(f.tm - r.tm);
    if (tmDiff > 2.5) continue;
    // Strategy-specific pair validation
    if (exons) {
      if (strategy === "junction") {
        const fSpan = spansJunction(f.start, f.end);
        const rSpan = spansJunction(r.bindStart, r.bindEnd);
        if (!fSpan && !rSpan) continue;
      }
      if (strategy === "neighboring") {
        if (selectedExons.length === 2) {
          const [ex1, ex2] = selectedExons;
          // Forward primer must be primarily on ex1, reverse on ex2
          if (f.exonNum !== ex1 || r.exonNum !== ex2) continue;
        } else {
          // No specific selection: just require fwd and rev on different exons
          if (f.exonNum === null || r.exonNum === null) continue;
          if (f.exonNum === r.exonNum) continue;
          // Also ensure fwd exon comes BEFORE rev exon (transcript order)
          if (f.exonNum > r.exonNum) continue;
        }
      }
    }
    const cross = dimerDG(f.seq, r.seq).dg;
    const score = tmDiff + Math.abs((f.gc + r.gc) / 2 - 50) / 5 +
      Math.abs(amp - (aMin + aMax) / 2) / 50 + Math.max(0, -cross - 6) / 2 +
      f.issues.length + r.issues.length;
    pairs.push({ forward: f, reverse: r, amplicon: amp, tmDiff, cross, score });
  }
  pairs.sort((a, b) => a.score - b.score);
  if (!pairs.length) {
    const stratMsg = strategy === "junction"
      ? " No junction-spanning primers found — try 'Neighboring exons' strategy or fetch a region with more exons."
      : strategy === "neighboring" && selectedExons.length === 2
        ? ` No pairs found between E${selectedExons[0]} and E${selectedExons[1]}. Try selecting adjacent exons.`
        : "";
    return { error: `No primer pair met the criteria.${stratMsg} Try widening the amplicon range or Tm in Advanced settings.` };
  }
  const top6 = pairs.slice(0, 6).map((pr) => ({
    ...pr,
    map: {
      length: seq.length, seqStr: seq,
      features: [
        { key: "forward", start: pr.forward.start, end: pr.forward.end, label: "forward", color: "#0d9488", strand: "fwd" },
        { key: "reverse", start: pr.reverse.bindStart, end: pr.reverse.bindEnd, label: "reverse", color: "#2563c9", strand: "rev" },
      ],
      spans: [{ start: pr.forward.start, end: pr.reverse.bindEnd, label: "amplicon", size: pr.amplicon }],
    },
  }));
  return { method: "Standard qPCR / PCR", pairs: top6, p,
           exonStrategy: strategy, exonsUsed: exons ? true : false };
}

/* ============================================================
   TETRA-PRIMER ARMS-PCR  (SNP genotyping, Ye et al. 2001)
   Joint optimizer: extension-aware deliberate-mismatch selection
   (Huang 1992 / Simsek & Adnan 2000) + joint Tm balance + band
   geometry. Allele→inner assignment searched both ways to escape
   weak-strand traps (C·T leak, terminal-A weakness).
   ============================================================ */

// --- Empirical extension-rate model (Huang 1992) ---
// -log10 relative extension efficiency for primer.template mispairs.
// Higher = more suppressed = better discrimination for ARMS.
const HUANG = {
  "A.C":3.5,"C.A":3.5,"G.T":3.5,"T.G":3.5,
  "T.C":4.5,"T.T":4.5,"A.A":6.0,
  "A.G":6.5,"G.A":6.5,"G.G":6.5,"C.C":6.5,
  "C.T":2.0, // C.T is the easy-extension exception — warns in UI
};
const PENULT_ADDED = 1.0;
const DOUBLE_MM_BONUS = 0.5;
const PENULT_MATCHED_PENALTY = { A:1.3, T:1.0, C:0.6, G:0.6 };
const MATCHED_LOSS_BUDGET = 1.6;

function mispairSuppression(primerBase, templateBase) {
  if (templateBase === COMP[primerBase]) return 0.0;
  return HUANG[`${primerBase}.${templateBase}`] ?? 3.0;
}

function evalDiscrimination(detected, other, orientation, deliberateM2, genomicM2Sense) {
  // forward primer 3' base = detected; reverse = COMP[detected]
  const pTerm = orientation === "forward" ? detected : COMP[detected];
  const tMatch = orientation === "forward" ? COMP[detected] : detected;
  const tWrong = orientation === "forward" ? COMP[other] : other;
  const termSup = mispairSuppression(pTerm, tWrong);

  // weak-discriminator flags
  const warnings = [];
  if (pTerm === "A") warnings.push(`3′-terminal A is a weak discriminator (Simsek & Adnan); consider the other strand.`);
  if (`${pTerm}.${tWrong}` === "C.T") warnings.push(`Wrong-allele mispair is C·T — Taq extends this ~10⁻² (Huang 1992); poor discrimination. Consider the other strand.`);

  const templateM2 = orientation === "forward" ? COMP[genomicM2Sense] : genomicM2Sense;
  const deliberateIsMismatch = (templateM2 !== COMP[deliberateM2]);
  let added = 0.0, matchedLoss = 0.0;
  if (deliberateIsMismatch) {
    added += PENULT_ADDED;
    if (termSup > 0) added += DOUBLE_MM_BONUS;
    matchedLoss = PENULT_MATCHED_PENALTY[deliberateM2] ?? 0.8;
  } else {
    warnings.push(`Chosen −2 base matches the template; no additional ARMS discrimination from deliberate mismatch.`);
  }
  const discriminationLog = termSup + added;
  const primable = matchedLoss <= MATCHED_LOSS_BUDGET;
  return { discriminationLog, matchedLoss, primable, warnings, pTerm, tWrong };
}

function bestDeliberate(detected, other, orientation, genomicM2Sense) {
  const templateM2 = orientation === "forward" ? COMP[genomicM2Sense] : genomicM2Sense;
  let best = null;
  for (const d2 of ["A","C","G","T"]) {
    if (COMP[d2] === templateM2) continue; // must be a real mismatch at -2
    const ev = evalDiscrimination(detected, other, orientation, d2, genomicM2Sense);
    if (!ev.primable) continue;
    if (!best || ev.discriminationLog > best.ev.discriminationLog) best = { d2, ev };
  }
  if (!best) {
    // fallback: even non-primable, pick best discrimination
    for (const d2 of ["C","G","A","T"]) {
      if (COMP[d2] === templateM2) continue;
      const ev = evalDiscrimination(detected, other, orientation, d2, genomicM2Sense);
      if (!best || ev.discriminationLog > best.ev.discriminationLog) best = { d2, ev };
    }
  }
  return best ?? { d2: "C", ev: evalDiscrimination(detected, other, orientation, "C", genomicM2Sense) };
}

function buildFI(full, snp, detected, L) {
  if (snp - L + 1 < 0) return null;
  const genomicM2Sense = full[snp - 1];
  const { d2, ev } = bestDeliberate(detected, detected === "A" ? "G" : "A", "forward", genomicM2Sense);
  const body = full.substring(snp - L + 1, snp - 1);
  const seq = body + d2 + detected;
  return { seq, d2, ev, Lfi: L };
}

function buildRI(full, snp, detected, L) {
  if (snp + L > full.length) return null;
  const genomicM2Sense = full[snp + 1];
  const { d2, ev } = bestDeliberate(detected, detected === "A" ? "G" : "A", "reverse", genomicM2Sense);
  const window = full.substring(snp, snp + L);
  const rc = revComp(window);
  const seq = rc.slice(0, -2) + d2 + COMP[detected];
  return { seq, d2, ev, Lri: L };
}

function designTetraPrimer(raw, p) {
  const norm = (raw || "").toUpperCase().replace(/U/g, "T");
  const match = norm.match(/\[([ACGT])\/([ACGT])\]/);
  if (!match) return { error: "Mark the SNP in brackets with both alleles, e.g. …GCT[A/G]TCA…" };
  const left = cleanDNA(norm.slice(0, match.index));
  const right = cleanDNA(norm.slice(match.index + match[0].length));
  const a1 = match[1], a2 = match[2];
  if (left.length < 80 || right.length < 80)
    return { error: "Provide ≥80 nt of flanking sequence each side of the SNP (≥150 nt recommended so the two allele bands separate well)." };
  const full = left + a1 + right;
  const snp = left.length;

  // outer primer search (unchanged geometry), memoized: the search only
  // depends on (Lfi, anchor) for fwd and (Lri, anchor) for rev, so cache by
  // those keys — otherwise it re-runs ~400 evalPrimer calls inside every
  // allele x Lfi x Lri x size x orientation combination (millions of calls,
  // which froze the UI on long sequences).
  const _fwdCache = new Map();
  const _revCache = new Map();
  const bestFwd = (snpPos, Lfi, anchor) => {
    const key = Lfi + ":" + anchor;
    if (_fwdCache.has(key)) return _fwdCache.get(key);
    let best = null;
    for (let s = Math.max(0, anchor - 25); s <= anchor + 25; s++)
      for (let L = 18; L <= 25; L++) {
        if (s + L > snpPos - Lfi) continue;
        const e = evalPrimer(full.substr(s, L), p);
        if (e.tm == null) continue;
        const sc = Math.abs(e.tm - p.target) + e.issues.length * 0.5;
        if (!best || sc < best.sc) best = { ...e, start: s, sc };
      }
    _fwdCache.set(key, best);
    return best;
  };
  const bestRev = (snpPos, Lri, anchor) => {
    const key = Lri + ":" + anchor;
    if (_revCache.has(key)) return _revCache.get(key);
    let best = null;
    for (let end = Math.min(full.length, anchor + 25); end >= anchor - 25; end--)
      for (let L = 18; L <= 25; L++) {
        const start = end - L;
        if (start < snpPos + Lri || end > full.length) continue;
        const e = evalPrimer(revComp(full.substr(start, L)), p);
        if (e.tm == null) continue;
        const sc = Math.abs(e.tm - p.target) + e.issues.length * 0.5;
        if (!best || sc < best.sc) best = { ...e, bindEnd: end, sc };
      }
    _revCache.set(key, best);
    return best;
  };

  // joint search: vary inner lengths + allele assignment, score all four Tm together
  let bestDesign = null, bestScore = Infinity;
  const SIZE_PAIRS = [];
  for (let sm = 100; sm <= 280; sm += 30)
    for (let bg = sm + 50; bg <= 310; bg += 30)
      SIZE_PAIRS.push([sm, bg]);

  for (const [fiAllele, riAllele] of [[a1, a2], [a2, a1]]) {
    for (let Lfi = 18; Lfi <= 27; Lfi += 2) {
      const fiR = buildFI(full, snp, fiAllele, Lfi);
      if (!fiR) continue;
      for (let Lri = 18; Lri <= 27; Lri += 2) {
        const riR = buildRI(full, snp, riAllele, Lri);
        if (!riR) continue;
        for (const [sA, sB] of SIZE_PAIRS) {
          for (const [sizeRI, sizeFI] of [[sA, sB], [sB, sA]]) {
            const a = snp - (sizeRI - Lri);
            const b = snp + (sizeFI - Lfi);
            if (a < 0 || b >= full.length) continue;
            const FO = bestFwd(snp, Lfi, a);
            const RO = bestRev(snp, Lri, b);
            if (!FO || !RO) continue;
            const tms = [FO.tm, tm(fiR.seq, p), tm(riR.seq, p), RO.tm].filter(Boolean);
            if (tms.length < 4) continue;
            const tmSpread = Math.max(...tms) - Math.min(...tms);
            const minDisc = Math.min(fiR.ev.discriminationLog, riR.ev.discriminationLog);
            const bandSep = Math.abs((snp + Lri - FO.start) - (RO.bindEnd - (snp - Lfi + 1)));
            const score = 1.4 * tmSpread + 1.0 * Math.max(0, 6 - minDisc) + 1.2 * Math.max(0, 40 - bandSep);
            if (score < bestScore) {
              bestScore = score;
              bestDesign = { FO, RO, fiR, riR, Lfi, Lri, fiAllele, riAllele };
            }
          }
        }
      }
    }
  }

  if (!bestDesign) return { error: "Could not place all four primers. Provide more flanking sequence or adjust the target Tm." };

  const { FO, RO, fiR, riR, Lfi, Lri, fiAllele, riAllele } = bestDesign;
  const fi = fiR.seq, ri = riR.seq;
  const len1 = snp + Lri - FO.start;       // FO + RI → riAllele amplicon
  const len2 = RO.bindEnd - (snp - Lfi + 1); // FI + RO → fiAllele amplicon
  const control = RO.bindEnd - FO.start;

  // collect warnings
  const allWarnings = [...new Set([...fiR.ev.warnings, ...riR.ev.warnings])];
  if (Math.abs(len1 - len2) < 40) allWarnings.push("The two allele bands differ by under 40 bp — add more flanking sequence for cleaner gel resolution.");
  const warn = allWarnings.length ? allWarnings.join(" ") : null;

  // discrimination summary for display
  const discFI = fiR.ev.discriminationLog.toFixed(1);
  const discRI = riR.ev.discriminationLog.toFixed(1);

  return {
    method: "Tetra-primer ARMS-PCR (SNP genotyping)",
    blurb: `Two allele-specific inner primers with jointly optimised Tm balance and extension-aware deliberate −2 mismatch selection (Huang 1992 / Simsek & Adnan 2000), plus two outer control primers (Ye et al. 2001). Each allele yields a distinct band; outer pair confirms amplification. Predicted wrong-allele suppression: FI ${discFI} log-fold · RI ${discRI} log-fold.`,
    arms: { a1, a2 },
    oligos: [
      { label: `Forward inner — allele ${fiAllele}`, seq: fi,
        segs: [[fi.slice(0,-2),"mir"],[fi.slice(-2,-1),"overhang"],[fi.slice(-1),"tail"]],
        evalIt: true,
        note: `3′ base = ${fiAllele} (allele-specific); deliberate −2 mismatch = ${fiR.d2}; predicted suppression ${discFI} log-fold (Huang 1992).` },
      { label: `Reverse inner — allele ${riAllele}`, seq: ri,
        segs: [[ri.slice(0,-2),"mir"],[ri.slice(-2,-1),"overhang"],[ri.slice(-1),"tail"]],
        evalIt: true,
        note: `3′ base = ${COMP[riAllele]} (complement of ${riAllele}); deliberate −2 mismatch = ${riR.d2}; predicted suppression ${discRI} log-fold (Huang 1992).` },
      { label: "Forward outer (control)", seq: FO.seq, segs: [[FO.seq,"scaffold"]], evalIt: true },
      { label: "Reverse outer (control)", seq: RO.seq, segs: [[RO.seq,"scaffold"]], evalIt: true },
    ],
    bands: [
      { name: `Allele ${riAllele}  —  forward outer + reverse inner`, size: len1 },
      { name: `Allele ${fiAllele}  —  forward inner + reverse outer`, size: len2 },
      { name: "Control  —  outer pair", size: control },
      { name: "Band separation", size: Math.abs(len1 - len2), note: "≥40 bp recommended" },
    ],
    map: {
      length: full.length, snp, seqStr: full,
      features: [
        { key: "Forward outer (control)", start: FO.start, end: FO.start + FO.seq.length, label: "F outer", color: "#64748b", strand: "fwd" },
        { key: `Forward inner — allele ${fiAllele}`, start: snp - Lfi + 1, end: snp + 1, label: `F inner · ${fiAllele}`, color: "#0d9488", strand: "fwd" },
        { key: `Reverse inner — allele ${riAllele}`, start: snp, end: snp + Lri, label: `R inner · ${riAllele}`, color: "#d97706", strand: "rev" },
        { key: "Reverse outer (control)", start: RO.bindEnd - RO.seq.length, end: RO.bindEnd, label: "R outer", color: "#64748b", strand: "rev" },
      ],
      spans: [
        { start: FO.start, end: snp + Lri, label: `allele ${riAllele}`, size: len1 },
        { start: snp - Lfi + 1, end: RO.bindEnd, label: `allele ${fiAllele}`, size: len2 },
        { start: FO.start, end: RO.bindEnd, label: "control", size: control },
      ],
    },
    warn,
    p,
  };
}

/* ============================================================
   RESTRICTION SITES
   A compact panel of common, widely-stocked enzymes. Recognition
   sequences may include standard IUPAC ambiguity codes (R/Y/N
   etc.); both strands are scanned since most sites are palindromic
   or near-palindromic.
   ============================================================ */

const ENZYMES = [
  { name: "EcoRI", site: "GAATTC" },
  { name: "BamHI", site: "GGATCC" },
  { name: "HindIII", site: "AAGCTT" },
  { name: "NotI", site: "GCGGCCGC" },
  { name: "XhoI", site: "CTCGAG" },
  { name: "PstI", site: "CTGCAG" },
  { name: "SalI", site: "GTCGAC" },
  { name: "NcoI", site: "CCATGG" },
  { name: "KpnI", site: "GGTACC" },
  { name: "SacI", site: "GAGCTC" },
  { name: "HaeIII", site: "GGCC" },
  { name: "MspI", site: "CCGG" },
  { name: "AluI", site: "AGCT" },
  { name: "TaqI", site: "TCGA" },
  { name: "RsaI", site: "GTAC" },
];

const IUPAC = { A: "A", C: "C", G: "G", T: "T", R: "[AG]", Y: "[CT]", N: "[ACGT]", W: "[AT]", S: "[GC]", K: "[GT]", M: "[AC]" };
const siteToRegex = (site) => new RegExp(site.split("").map((b) => IUPAC[b] || b).join(""), "g");

function findRestrictionSites(seqStr, enzymeNames) {
  if (!seqStr || !enzymeNames?.length) return [];
  const hits = [];
  for (const enz of ENZYMES) {
    if (!enzymeNames.includes(enz.name)) continue;
    const re = siteToRegex(enz.site);
    let m;
    while ((m = re.exec(seqStr))) hits.push({ enzyme: enz.name, start: m.index, end: m.index + enz.site.length });
    const rc = revComp(enz.site);
    if (rc !== enz.site) {
      const reR = siteToRegex(rc);
      while ((m = reR.exec(seqStr))) hits.push({ enzyme: enz.name, start: m.index, end: m.index + rc.length });
    }
  }
  return hits.sort((a, b) => a.start - b.start);
}

/* ============================================================
   EXPORT
   Normalizes any of the three result shapes into a flat list of
   real, orderable oligos (skips placeholder text like the kit
   adapter line in poly(A) mode), then offers CSV / FASTA / SVG.
   ============================================================ */

function getExportOligos(result, selectedPair = 0) {
  if (!result) return [];
  // miRNA / standard design: oligos array with segs
  if (result.oligos) {
    return result.oligos
      .filter((o) => o.segs && o.segs.length > 0)
      .map((o) => {
        const seq = cleanDNA(o.seq);
        const e = o.evalIt ? evalPrimer(seq, result.p) : null;
        return { name: o.label, seq, tm: e?.tm ?? null, gc: e?.gc ?? null, length: seq.length };
      });
  }
  // Standard PCR / SNP-aware: pairs array
  if (result.pairs) {
    const pr = result.pairs[selectedPair] || result.pairs[0];
    return [
      { name: "Forward primer", seq: pr.forward.seq, tm: pr.forward.tm, gc: pr.forward.gc, length: pr.forward.seq.length },
      { name: "Reverse primer", seq: pr.reverse.seq, tm: pr.reverse.tm, gc: pr.reverse.gc, length: pr.reverse.seq.length },
    ];
  }
  // HRM: primers array (each has .seq .role .tm .gc)
  if (result.primers && Array.isArray(result.primers) && result.primers[0]?.role !== undefined) {
    return result.primers
      .filter(p => p.seq && p.seq.length > 0)
      .map(p => ({
        name: p.role ? `${p.role.charAt(0).toUpperCase() + p.role.slice(1)} primer${p.note ? ` (${p.note.slice(0,30)})` : ""}` : "Primer",
        seq: p.seq, tm: p.tm ?? null, gc: p.gc ?? null, length: p.seq.length,
      }));
  }
  // Methylation: sets object {msp: [...], bsp: [...], ms_hrm: [...]}
  if (result.sets && typeof result.sets === "object") {
    return Object.entries(result.sets).flatMap(([setName, primers]) =>
      (primers || []).filter(p => p.seq).map(p => ({
        name: `${setName.toUpperCase()} — ${p.role || "primer"}`,
        seq: p.seq, tm: p.tm ?? null, gc: p.gc ?? null, length: p.seq.length,
      }))
    );
  }
  return [];
}

function download(filename, content, mime) {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

function exportCSV(oligos, methodName) {
  const rows = [["Name", "Sequence (5'->3')", "Length (nt)", "Tm (C)", "GC (%)"]];
  oligos.forEach((o) => rows.push([o.name, o.seq, o.length, o.tm != null ? o.tm.toFixed(1) : "", o.gc != null ? o.gc.toFixed(0) : ""]));
  const csv = rows.map((r) => r.map((c) => `"${String(c).replace(/"/g, '""')}"`).join(",")).join("\r\n");
  download(`${(methodName || "primers").replace(/[^\w-]+/g, "_")}.csv`, csv, "text/csv");
}

function exportFASTA(oligos, methodName) {
  const fasta = oligos.map((o) => `>${o.name.replace(/\s+/g, "_")}\n${o.seq}`).join("\n");
  download(`${(methodName || "primers").replace(/[^\w-]+/g, "_")}.fasta`, fasta, "text/plain");
}

function exportSVG(svgEl, filename) {
  if (!svgEl) return;
  const clone = svgEl.cloneNode(true);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  const bg = document.createElementNS("http://www.w3.org/2000/svg", "rect");
  bg.setAttribute("x", "0"); bg.setAttribute("y", "0"); bg.setAttribute("width", "100%"); bg.setAttribute("height", "100%"); bg.setAttribute("fill", "#ffffff");
  clone.insertBefore(bg, clone.firstChild);
  const src = new XMLSerializer().serializeToString(clone);
  download(filename, src, "image/svg+xml");
}

function ExportBar({ result, mapId, selectedPair = 0, onSave, justSaved }) {
  const oligos = useMemo(() => getExportOligos(result, selectedPair), [result, selectedPair]);
  if (!oligos.length) return null;
  const name = (result.method || "primers").toLowerCase().replace(/[^\w-]+/g, "_");
  return (
    <div className="mt-4 flex flex-wrap items-center gap-2 rounded-lg border border-slate-200 bg-slate-50 p-3">
      <span className="text-xs font-medium text-slate-500 mr-1">Export</span>
      <button onClick={() => exportCSV(oligos, name)}
        className="inline-flex items-center gap-1.5 rounded-md bg-white px-2.5 py-1.5 text-xs font-medium text-slate-700 ring-1 ring-inset ring-slate-200 hover:bg-slate-100">
        <Download size={12} /> Order sheet (CSV)
      </button>
      <button onClick={() => exportFASTA(oligos, name)}
        className="inline-flex items-center gap-1.5 rounded-md bg-white px-2.5 py-1.5 text-xs font-medium text-slate-700 ring-1 ring-inset ring-slate-200 hover:bg-slate-100">
        <Download size={12} /> FASTA
      </button>
      {mapId && (
        <button onClick={() => exportSVG(document.getElementById(mapId), `${name}_map.svg`)}
          className="inline-flex items-center gap-1.5 rounded-md bg-white px-2.5 py-1.5 text-xs font-medium text-slate-700 ring-1 ring-inset ring-slate-200 hover:bg-slate-100">
          <Download size={12} /> Map (SVG)
        </button>
      )}
      {onSave && (
        <button onClick={onSave}
          className={`ml-auto inline-flex items-center gap-1.5 rounded-md px-2.5 py-1.5 text-xs font-semibold transition-colors ${justSaved ? "bg-teal-600 text-white" : "bg-slate-900 text-white hover:bg-slate-800"}`}>
          {justSaved ? <><Check size={12} /> Saved</> : <><Save size={12} /> Save design</>}
        </button>
      )}
    </div>
  );
}

/* ============================================================
   SAVED DESIGNS WORKSPACE  (localStorage persistence)
   ============================================================ */

const STORAGE_KEY = "freeprimers.designs.v1";

function loadDesigns() {
  try { const raw = localStorage.getItem(STORAGE_KEY); return raw ? JSON.parse(raw) : []; }
  catch { return []; }
}
function persistDesigns(designs) {
  try { localStorage.setItem(STORAGE_KEY, JSON.stringify(designs)); } catch { /* quota/disabled */ }
}

function useSavedDesigns() {
  const [designs, setDesigns] = useState(() => loadDesigns());
  useEffect(() => { persistDesigns(designs); }, [designs]);
  const save = (entry) => {
    setDesigns((prev) => {
      const e = { ...entry, id: entry.id || `d_${Date.now()}_${Math.random().toString(36).slice(2, 7)}`, savedAt: Date.now() };
      const without = prev.filter((d) => d.id !== e.id);
      return [e, ...without].slice(0, 200);
    });
  };
  const remove = (id) => setDesigns((prev) => prev.filter((d) => d.id !== id));
  const clear = () => setDesigns([]);
  return { designs, save, remove, clear };
}

function timeAgo(ts) {
  const s = Math.floor((Date.now() - ts) / 1000);
  if (s < 60) return "just now";
  if (s < 3600) return `${Math.floor(s / 60)}m ago`;
  if (s < 86400) return `${Math.floor(s / 3600)}h ago`;
  const d = Math.floor(s / 86400);
  if (d < 30) return `${d}d ago`;
  return new Date(ts).toLocaleDateString();
}

const TAB_LABELS = {
  mirna: "miRNA", standard: "Standard qPCR", snp: "SNP genotyping",
  placement: "SNP-aware", hrm: "HRM", methyl: "Methylation",
};

function WorkspaceDrawer({ open, onClose, designs, onLoad, onRemove, onClear, onExportAll }) {
  if (!open) return null;
  return (
    <div className="fixed inset-0 z-50 flex justify-end" onClick={onClose}>
      <div className="absolute inset-0 bg-slate-900/20 backdrop-blur-[1px]" />
      <div onClick={(e) => e.stopPropagation()}
        className="relative h-full w-full max-w-md overflow-y-auto border-l border-slate-200 bg-white shadow-2xl">
        <div className="sticky top-0 z-10 flex items-center justify-between border-b border-slate-200 bg-white px-4 py-3">
          <div className="flex items-center gap-2">
            <Bookmark size={16} className="text-teal-600" />
            <span className="text-sm font-semibold text-slate-800">Saved designs</span>
            <span className="rounded-full bg-slate-100 px-2 py-0.5 text-[11px] font-medium text-slate-500">{designs.length}</span>
          </div>
          <button onClick={onClose} className="rounded-md p-1 text-slate-400 hover:bg-slate-100 hover:text-slate-700"><X size={16} /></button>
        </div>
        {designs.length > 0 && (
          <div className="flex items-center gap-2 border-b border-slate-100 px-4 py-2">
            <button onClick={onExportAll}
              className="inline-flex items-center gap-1.5 rounded-md bg-slate-900 px-2.5 py-1 text-[11px] font-semibold text-white hover:bg-slate-800">
              <Download size={11} /> Export panel (all oligos)
            </button>
            <button onClick={() => { if (confirm("Remove all saved designs? This cannot be undone.")) onClear(); }}
              className="ml-auto inline-flex items-center gap-1 text-[11px] text-slate-400 hover:text-rose-500">
              <Trash2 size={11} /> Clear all
            </button>
          </div>
        )}
        {designs.length === 0 ? (
          <div className="flex flex-col items-center justify-center px-6 py-16 text-center">
            <FolderOpen size={32} className="mb-3 text-slate-300" />
            <p className="text-sm font-medium text-slate-600">No saved designs yet</p>
            <p className="mt-1 text-[12px] text-slate-400">Design primers, then click <span className="font-medium text-slate-600">Save design</span> to keep them across sessions and build an order panel.</p>
          </div>
        ) : (
          <div className="divide-y divide-slate-100">
            {designs.map((d) => (
              <div key={d.id} className="group px-4 py-3 hover:bg-slate-50">
                <div className="flex items-start justify-between gap-2">
                  <button onClick={() => onLoad(d)} className="min-w-0 flex-1 text-left">
                    <div className="flex items-center gap-2">
                      <span className="truncate text-sm font-semibold text-slate-800">{d.name}</span>
                      <span className="shrink-0 rounded bg-teal-50 px-1.5 py-0.5 text-[10px] font-medium text-teal-700">{TAB_LABELS[d.tab] || d.tab}</span>
                    </div>
                    <div className="mt-0.5 flex flex-wrap items-center gap-x-2 gap-y-0.5 text-[11px] text-slate-400">
                      <span className="inline-flex items-center gap-1"><Clock size={10} /> {timeAgo(d.savedAt)}</span>
                      {d.oligoCount != null && <span>· {d.oligoCount} oligo{d.oligoCount !== 1 ? "s" : ""}</span>}
                      {d.summary && <span className="truncate">· {d.summary}</span>}
                    </div>
                  </button>
                  <div className="flex shrink-0 items-center gap-1 opacity-0 transition-opacity group-hover:opacity-100">
                    <button onClick={() => onLoad(d)} title="Load design" className="rounded p-1 text-slate-400 hover:bg-teal-50 hover:text-teal-600"><FolderOpen size={13} /></button>
                    <button onClick={() => onRemove(d.id)} title="Delete" className="rounded p-1 text-slate-400 hover:bg-rose-50 hover:text-rose-500"><Trash2 size={13} /></button>
                  </div>
                </div>
              </div>
            ))}
          </div>
        )}
        <p className="px-4 py-3 text-[10px] text-slate-300">Saved locally in this browser only — not uploaded anywhere. Clearing browser data removes them.</p>
      </div>
    </div>
  );
}

function exportPanelCSV(designs) {
  const rows = [["Design", "Type", "Oligo name", "Sequence (5'->3')", "Length", "Tm", "GC%", "Saved"]];
  designs.forEach((d) => {
    (d.oligos || []).forEach((o) => {
      rows.push([d.name, TAB_LABELS[d.tab] || d.tab, o.name, o.seq, o.length ?? o.seq?.length ?? "",
                 o.tm != null ? (o.tm.toFixed ? o.tm.toFixed(1) : o.tm) : "",
                 o.gc != null ? (o.gc.toFixed ? o.gc.toFixed(0) : o.gc) : "",
                 new Date(d.savedAt).toLocaleDateString()]);
    });
  });
  const csv = rows.map((r) => r.map((c) => `"${String(c).replace(/"/g, '""')}"`).join(",")).join("\n");
  const blob = new Blob([csv], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = `primer_panel_${designs.length}_designs.csv`;
  a.click(); URL.revokeObjectURL(url);
}

/* ============================================================
   UI
   ============================================================ */

const PRESETS = [
  { name: "hsa-miR-21-5p", seq: "UAGCUUAUCAGACUGAUGUUGA" },
  { name: "hsa-let-7a-5p", seq: "UGAGGUAGUAGGUUGUAUAGUU" },
  { name: "hsa-miR-16-5p", seq: "UAGCAGCACGUAAAUAUUGGCG" },
  { name: "hsa-miR-122-5p", seq: "UGGAGUGUGACAAUGGUGUUUG" },
];

const SNP_EXAMPLE =
  "ATGGCCTAGCATTACGGATCCAATGCCTGAGTCAGGTTACAGCCTAGGATCCGATTACGGCTAACGTTGACCTAGGCATTACGGATCCAATGCCTGAGTCAGGTTACAGCCTAGGATCCGATTACGGCTAACGTTGACCTAGGCATTACGGATCCAATGCCTGAGTCAGGTTACAGCCTAGGATCCGATTACGGCTAACGTTGACCTAGGCATTACGGATCCAATGCCTGAGTCAGGTTACAGCCTAGGATCCGATTAC" +
  "[C/T]" +
  "TGGATCAATGCCTAGGTTACAGCCTAGGATCCGATTACGGCTAACGTTGACCTAGGCATTACGGATCCAATGCCTAGGTTACAGCCTAGGATCCGATTACGGCTAACGTTGACCTAGGCATTACGGATCCAATGCCTAGGTTACAGCCTAGGATCCGATTACGGCTAACGTTGACCTAGGCATTACGGATCCAATGCCTAGGTTACAGCCTAGGATCCGATTACGGCTAACGTTGACCTAGGCATTACGGATCCAAT";

const SEG = {
  scaffold: "bg-slate-200 text-slate-700",
  overhang: "bg-rose-100 text-rose-700",
  tail: "bg-amber-100 text-amber-800",
  mir: "bg-teal-100 text-teal-800",
};
const SEG_LABEL = { scaffold: "universal scaffold", overhang: "3′ overhang", tail: "Tm-adjusting tail", mir: "miRNA-derived" };

const Q = {
  good: { c: "text-emerald-700 bg-emerald-50 ring-emerald-600/20", Icon: CircleCheck, t: "Good" },
  ok: { c: "text-amber-700 bg-amber-50 ring-amber-600/20", Icon: CircleAlert, t: "Acceptable" },
  warn: { c: "text-rose-700 bg-rose-50 ring-rose-600/20", Icon: TriangleAlert, t: "Review" },
};

function CopyBtn({ text }) {
  const [done, setDone] = useState(false);
  return (
    <button
      onClick={async () => {
        try { await navigator.clipboard.writeText(text); } catch { /* ignore */ }
        setDone(true); setTimeout(() => setDone(false), 1200);
      }}
      className="inline-flex items-center gap-1 rounded-md px-2 py-1 text-xs font-medium text-slate-500 hover:bg-slate-100 hover:text-slate-800 transition-colors"
    >
      {done ? <Check size={13} /> : <Copy size={13} />}{done ? "Copied" : "Copy"}
    </button>
  );
}

function Seq({ segs, seq }) {
  if (!segs || !segs.length) return <span className="text-slate-400 italic">{seq}</span>;
  return (
    <span className="break-all">
      {segs.map(([s, kind], i) => s ? <span key={i} className={`${SEG[kind]} rounded px-0.5 py-px`}>{s}</span> : null)}
    </span>
  );
}

function Metric({ label, value, hint }) {
  return (
    <div className="flex flex-col" title={hint}>
      <span className="text-[10px] uppercase tracking-wider text-slate-400">{label}</span>
      <span className="font-mono text-sm text-slate-800">{value}</span>
    </div>
  );
}

function OligoCard({ o, p, hoverKey, onHover, backendActive }) {
  const e = o.evalIt ? evalPrimer(cleanDNA(o.seq), p) : null;
  const q = e ? Q[e.quality] : null;
  const active = hoverKey === o.label;
  const [verified, setVerified] = useState(null);
  const [verifying, setVerifying] = useState(false);

  useEffect(() => {
    setVerified(null);
    if (!backendActive || !e) return;
    let live = true;
    setVerifying(true);
    backendVerify(cleanDNA(o.seq), p).then((res) => {
      if (live) { setVerified(res); setVerifying(false); }
    });
    return () => { live = false; };
  }, [backendActive, o.seq, p.na, p.mg, p.dntp, p.primer]);

  const shown = {
    tm: verified?.tm ?? e?.tm,
    hairpin: verified?.hairpin ?? e?.hairpin,
    selfDimer: verified?.selfDimer ?? e?.selfDimer,
  };

  return (
    <div onMouseEnter={() => onHover?.(o.label)} onMouseLeave={() => onHover?.(null)}
      className={`rounded-lg border p-4 transition-colors ${active ? "border-slate-400 bg-slate-50" : "border-slate-200 bg-white"}`}>
      <div className="flex items-start justify-between gap-3">
        <div className="min-w-0">
          <div className="flex items-center gap-2 flex-wrap">
            <span className="text-sm font-semibold text-slate-800">{o.label}</span>
            {q && (
              <span className={`inline-flex items-center gap-1 rounded-full px-2 py-0.5 text-[11px] font-medium ring-1 ring-inset ${q.c}`}>
                <q.Icon size={11} /> {q.t}
              </span>
            )}
            {verifying && (
              <span className="inline-flex items-center gap-1 text-[11px] font-medium text-slate-400">
                <Loader2 size={11} className="animate-spin" /> verifying
              </span>
            )}
            {verified && !verifying && (
              <span className="inline-flex items-center gap-1 rounded-full bg-teal-50 px-2 py-0.5 text-[11px] font-medium text-teal-700 ring-1 ring-inset ring-teal-600/20">
                <ShieldCheck size={11} /> backend-verified
              </span>
            )}
          </div>
          <div className="mt-1.5 font-mono text-sm leading-relaxed">
            5′– <Seq segs={o.segs} seq={o.seq} /> –3′
          </div>
        </div>
        {cleanDNA(o.seq).length > 0 && o.segs.length > 0 && <CopyBtn text={cleanDNA(o.seq)} />}
      </div>

      {e && (
        <div className="mt-3 grid grid-cols-3 gap-3 border-t border-slate-100 pt-3 sm:grid-cols-5">
          <Metric label="Length" value={`${cleanDNA(o.seq).length} nt`} />
          <Metric label="Tm" value={shown.tm != null ? `${shown.tm.toFixed(1)}°` : "—"} hint={verified?.tm != null ? "primer3-py nearest-neighbor (backend)" : "SantaLucia 1998 nearest-neighbor (estimate)"} />
          <Metric label="GC" value={`${e.gc.toFixed(0)}%`} />
          <Metric label="Self-dimer" value={`${shown.selfDimer != null ? shown.selfDimer.toFixed(1) : "—"}`} hint={verified?.selfDimer != null ? "primer3-py heterodimer ΔG, kcal/mol (backend)" : "Estimated ΔG (kcal/mol)"} />
          <Metric label="Hairpin" value={`${shown.hairpin != null ? shown.hairpin.toFixed(1) : "—"}`} hint={verified?.hairpin != null ? "primer3-py hairpin ΔG, kcal/mol (backend)" : "Estimated ΔG (kcal/mol)"} />
        </div>
      )}
      {e && e.issues.length > 0 && (
        <ul className="mt-2 space-y-0.5">
          {e.issues.map((it, i) => <li key={i} className="text-xs text-amber-700">• {it}</li>)}
        </ul>
      )}
      {o.note && <p className="mt-2 text-xs text-slate-500">{o.note}</p>}
    </div>
  );
}

// Shared exon-membership computation so PairCard and PrimerAlignment never disagree.
// Returns { label, exonNums, spansJunction } for a primer footprint [pStart, pEnd).
// In cDNA mode exons are contiguous, so spanning two exons = a real splice junction.
function exonMembership(exons, pStart, pEnd) {
  if (!exons || pStart == null || pEnd == null) return null;
  const hits = exons.filter(ex => pStart < ex.end && pEnd > ex.start);
  if (!hits.length) return { label: "intron ⚠", exonNums: [], spansJunction: false };
  if (hits.length === 1) return { label: `E${hits[0].exon_number}`, exonNums: [hits[0].exon_number], spansJunction: false };
  // Spans 2+ exons. List them in order; this is a junction primer (good for cDNA).
  const nums = hits.map(h => h.exon_number).sort((a, b) => a - b);
  return { label: `E${nums.join("/E")} junction`, exonNums: nums, spansJunction: true };
}

function PairCard({ pair, idx, active, onSelect, hoverKey, onHover, backendActive, salts, exons }) {
  const { forward: f, reverse: r } = pair;
  const [crossDimer, setCrossDimer] = useState(null);
  const [checking, setChecking] = useState(false);

  // Use the shared helper so labels match the alignment exactly
  const fMem = exons ? exonMembership(exons, f.start, f.end) : null;
  const rMem = exons ? exonMembership(exons, r.bindStart, r.bindEnd) : null;
  const fExons = fMem?.label ?? null;
  const rExons = rMem?.label ?? null;
  // A pair "spans a junction" if either primer itself spans one, OR the two
  // primers sit on different exons (neighboring-exon design).
  const fNum = fMem?.exonNums?.[0];
  const rNum = rMem?.exonNums?.[0];
  const spansJunction = (fMem?.spansJunction || rMem?.spansJunction) ||
    (fMem && rMem && !fMem.label.includes("intron") && !rMem.label.includes("intron") && fNum !== rNum);

  useEffect(() => {
    setCrossDimer(null);
    if (!backendActive || !active) return;
    let live = true;
    setChecking(true);
    backendDuplex(f.seq, r.seq, salts).then((dg) => {
      if (live) { setCrossDimer(dg); setChecking(false); }
    });
    return () => { live = false; };
  }, [backendActive, active, f.seq, r.seq, salts?.na, salts?.mg, salts?.dntp, salts?.primer]);

  return (
    <div role="button" tabIndex={0} onClick={onSelect} onKeyDown={(e) => (e.key === "Enter" || e.key === " ") && onSelect()}
      className={`w-full cursor-pointer rounded-lg border p-4 text-left transition-colors ${active ? "border-teal-500 bg-teal-50/40 ring-1 ring-teal-500/30" : "border-slate-200 bg-white hover:border-slate-300"}`}>
      <div className="mb-3 flex items-center justify-between">
        <span className={`text-sm font-semibold ${active ? "text-teal-800" : "text-slate-800"}`}>
          Pair {idx + 1}{active && <span className="ml-2 text-[11px] font-medium text-teal-600">shown on map</span>}
        </span>
        <div className="flex items-center gap-4 text-xs text-slate-500">
          <span>Amplicon <span className="font-mono text-slate-800">{pair.amplicon} bp</span></span>
          <span>ΔTm <span className="font-mono text-slate-800">{pair.tmDiff.toFixed(1)}°</span></span>
          {spansJunction && <span className="rounded-full bg-teal-100 px-2 py-0.5 text-[10px] font-semibold text-teal-700">Junction-spanning ✓</span>}
        </div>
      </div>
      {[["Forward", f, "forward"], ["Reverse", r, "reverse"]].map(([lab, x, key]) => (
        <div key={lab} className={`mb-2 rounded-md px-1.5 py-1 -mx-1.5 last:mb-0 ${active && hoverKey === key ? "bg-slate-200/60" : ""}`}
          onMouseEnter={() => active && onHover?.(key)} onMouseLeave={() => active && onHover?.(null)}>
          <div className="flex items-center justify-between">
            <span className="text-xs font-medium text-slate-500">{lab}</span>
            <span onClick={(e) => e.stopPropagation()}><CopyBtn text={x.seq} /></span>
          </div>
          <div className="font-mono text-sm text-slate-800 break-all">5′–{x.seq}–3′</div>
          <div className="mt-0.5 flex flex-wrap gap-x-4 gap-y-0.5 text-[11px] text-slate-500">
            <span>Tm {x.tm.toFixed(1)}°</span><span>GC {x.gc.toFixed(0)}%</span><span>{x.seq.length} nt</span>
            {exons && (key === "forward" ? fExons : rExons) && (() => {
              const lbl = key === "forward" ? fExons : rExons;
              const isIntron = lbl.includes("intron");
              const isJunction = lbl.includes("junction");
              return (
                <span className={`font-medium ${isIntron ? "text-rose-600" : isJunction ? "text-teal-700" : "text-teal-600"}`}>
                  {isIntron ? "⚠ intron (gDNA only)" : lbl}
                </span>
              );
            })()}
          </div>
        </div>
      ))}
      {active && backendActive && (
        <div className="mt-2 flex items-center gap-1.5 border-t border-teal-200/60 pt-2 text-[11px]">
          {checking ? (
            <span className="inline-flex items-center gap-1 text-slate-400"><Loader2 size={11} className="animate-spin" /> checking forward×reverse cross-dimer…</span>
          ) : crossDimer != null ? (
            <span className={`inline-flex items-center gap-1 font-medium ${crossDimer < -6 ? "text-rose-600" : "text-teal-700"}`}>
              {crossDimer < -6 ? <ShieldAlert size={11} /> : <ShieldCheck size={11} />}
              Cross-dimer (F×R) ΔG {crossDimer.toFixed(1)} kcal/mol (backend)
            </span>
          ) : null}
        </div>
      )}
    </div>
  );
}

function AnnotationMap({ map, id, hoverKey, onHoverFeature, sites = [] }) {
  const { length, features = [], spans = [], snp = null, seqStr = null } = map;
  const W = 700, padX = 52, usable = W - padX * 2;
  const xOf = (pos) => padX + (length ? (pos / length) * usable : 0);
  const showBases = seqStr && length <= 60;
  const baseY = 92, fwdTop = 58, revTop = 100, blockH = 16;
  const sitesTop = revTop + blockH + (sites.length ? 26 : 4);
  const spansTop = sitesTop + (sites.length ? 22 : 0);
  const H = spansTop + spans.length * 26 + 14;
  const SANS = "ui-sans-serif, system-ui, sans-serif";

  const feat = (f, i) => {
    const active = hoverKey && f.key === hoverKey;
    const x0 = xOf(f.start), w = Math.max(4, xOf(f.end) - x0);
    const top = f.strand === "rev" ? revTop : fwdTop, ah = 6;
    const head = f.strand === "fwd"
      ? `${x0 + w},${top - 2} ${x0 + w + ah},${top + blockH / 2} ${x0 + w},${top + blockH + 2}`
      : `${x0},${top - 2} ${x0 - ah},${top + blockH / 2} ${x0},${top + blockH + 2}`;
    return (
      <g key={i} style={{ cursor: f.key ? "pointer" : "default" }}
        onMouseEnter={() => f.key && onHoverFeature?.(f.key)} onMouseLeave={() => f.key && onHoverFeature?.(null)}>
        {active && <rect x={x0 - 3} y={top - 3} width={w + 6} height={blockH + 6} rx="4" fill="none" stroke="#0f172a" strokeWidth="1.5" />}
        <rect x={x0} y={top} width={w} height={blockH} rx="2" fill={f.color} opacity={hoverKey && !active ? 0.45 : 0.92} />
        {f.strand && <polygon points={head} fill={f.color} opacity={hoverKey && !active ? 0.45 : 1} />}
        <text x={x0 + w / 2} y={f.strand === "rev" ? top + blockH + 12 : top - 4} fontSize="11" textAnchor="middle" fill={active ? "#0f172a" : "#475569"} fontWeight={active ? "600" : "400"} fontFamily={SANS}>{f.label}</text>
      </g>
    );
  };

  return (
    <svg id={id} viewBox={`0 0 ${W} ${H}`} width="100%" style={{ maxWidth: "100%", height: "auto" }} role="img" aria-label="primer map">
      <line x1={padX} y1={baseY} x2={W - padX} y2={baseY} stroke="#cbd5e1" strokeWidth="2" />
      {showBases
        ? seqStr.split("").map((b, i) => (
            <text key={`b${i}`} x={xOf(i + 0.5)} y={baseY + 4} fontSize={Math.min(13, usable / seqStr.length)} textAnchor="middle" fontFamily="ui-monospace, monospace" fill="#334155">{b}</text>
          ))
        : [0, 0.25, 0.5, 0.75, 1].map((fr, i) => (
            <g key={`t${i}`}>
              <line x1={xOf(length * fr)} y1={baseY - 4} x2={xOf(length * fr)} y2={baseY + 4} stroke="#94a3b8" />
              <text x={xOf(length * fr)} y={baseY + 18} fontSize="10" textAnchor="middle" fill="#94a3b8" fontFamily={SANS}>{Math.round(length * fr)}</text>
            </g>
          ))}
      {snp != null && (
        <g>
          <line x1={xOf(snp + 0.5)} y1={fwdTop - 8} x2={xOf(snp + 0.5)} y2={revTop + blockH + 8} stroke="#e11d48" strokeWidth="1.5" strokeDasharray="3 2" />
          <rect x={xOf(snp + 0.5) - 4} y={baseY - 4} width="8" height="8" fill="#e11d48" transform={`rotate(45 ${xOf(snp + 0.5)} ${baseY})`} />
          <text x={xOf(snp + 0.5)} y={fwdTop - 12} fontSize="10" textAnchor="middle" fill="#e11d48" fontFamily={SANS}>SNP</text>
        </g>
      )}
      {features.map(feat)}
      {sites.length > 0 && (
        <g>
          <line x1={padX} y1={sitesTop} x2={W - padX} y2={sitesTop} stroke="#e2e8f0" strokeWidth="1" />
          {sites.map((s, i) => {
            const x = xOf((s.start + s.end) / 2);
            return (
              <g key={`site${i}`}>
                <line x1={x} y1={sitesTop - 6} x2={x} y2={sitesTop + 6} stroke="#7c3aed" strokeWidth="1.5" />
                <text x={x} y={sitesTop + 18} fontSize="9.5" textAnchor="middle" fill="#7c3aed" fontFamily={SANS}>{s.enzyme}</text>
              </g>
            );
          })}
        </g>
      )}
      {spans.map((s, i) => {
        const y = spansTop + i * 26, x0 = xOf(s.start), x1 = xOf(s.end);
        return (
          <g key={`s${i}`}>
            <line x1={x0} y1={y} x2={x1} y2={y} stroke="#64748b" strokeWidth="1.5" />
            <line x1={x0} y1={y - 4} x2={x0} y2={y + 4} stroke="#64748b" strokeWidth="1.5" />
            <line x1={x1} y1={y - 4} x2={x1} y2={y + 4} stroke="#64748b" strokeWidth="1.5" />
            <text x={(x0 + x1) / 2} y={y - 6} fontSize="11" textAnchor="middle" fill="#334155" fontFamily={SANS}>{s.label}{s.size != null ? `  ·  ${s.size} bp` : ""}</text>
          </g>
        );
      })}
    </svg>
  );
}

// miRNA database search: calls /engines/mirna/fetch, shows results as cards,
// loads chosen sequence into the design textarea with one click.
function MirnaDbSearch({ onSelect }) {
  const [query, setQuery] = useState("");
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");
  const [selected, setSelected] = useState(null);
  const debounce = useRef(null);

  const ARM_COLOR = { "5p": "bg-teal-100 text-teal-700", "3p": "bg-violet-100 text-violet-700", "": "bg-slate-100 text-slate-500" };

  const search = async (q) => {
    if (!q.trim()) { setResults(null); return; }
    setLoading(true); setError("");
    try {
      const r = await fetch(`${getBackendUrl()}/engines/mirna/fetch`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ query: q.trim(), max_results: 12 }),
      });
      const j = await r.json();
      setResults(j.ok ? j.hits : []);
      if (!j.ok) setError(j.note || "No results.");
    } catch { setError("Backend not reachable."); setResults([]); }
    finally { setLoading(false); }
  };

  const onInput = (v) => {
    setQuery(v); setSelected(null);
    clearTimeout(debounce.current);
    debounce.current = setTimeout(() => search(v), 380);
  };

  const pick = (hit) => {
    setSelected(hit.name);
    onSelect(hit.sequence, hit);   // sequence + full hit for info card
  };

  return (
    <div className="mt-2 mb-1">
      {/* Search row — full width */}
      <div className="relative">
        <input
          value={query} onChange={(e) => onInput(e.target.value)}
          placeholder="Search miRBase 22.1 — e.g. let-7a, miR-21, MIMAT0000076…"
          className="w-full rounded-lg border border-slate-200 bg-white py-2 pl-3 pr-8 text-[12px] text-slate-700 outline-none focus:border-teal-500 focus:ring-2 focus:ring-teal-500/20"
        />
        {loading && <Loader2 size={14} className="absolute right-2.5 top-1/2 -translate-y-1/2 animate-spin text-teal-500" />}
        {!loading && query && (
          <button onClick={() => { setQuery(""); setResults(null); setError(""); }}
            className="absolute right-2.5 top-1/2 -translate-y-1/2 text-slate-400 hover:text-slate-600">
            <X size={14} />
          </button>
        )}
      </div>

      {/* Presets + FASTA row */}
      <div className="mt-2 flex flex-wrap items-center gap-2">
        <span className="text-[11px] text-slate-400">or try:</span>
        {PRESETS.map((pr) => (
          <button key={pr.name} onClick={() => { onSelect(pr.seq); setSelected(pr.name); setQuery(""); setResults(null); }}
            className="rounded-md bg-slate-100 px-2.5 py-1 font-mono text-[11px] text-slate-600 hover:bg-slate-200">
            {pr.name}
          </button>
        ))}
        <label className="ml-auto inline-flex cursor-pointer items-center gap-1 rounded-md bg-slate-100 px-2.5 py-1 text-[11px] text-slate-600 hover:bg-slate-200">
          <Upload size={12} /> Upload FASTA
          <input type="file" accept=".fa,.fasta,.txt" className="hidden" onChange={(e) => {
            const f = e.target.files?.[0]; if (!f) return;
            const rd = new FileReader();
            rd.onload = () => {
              const lines = String(rd.result || "").split(/\r?\n/);
              let s = ""; for (const ln of lines) { if (ln.startsWith(">")) { if (s) break; continue; } s += ln.trim(); }
              onSelect(s.toUpperCase()); setSelected("(uploaded)");
            };
            rd.readAsText(f); e.target.value = "";
          }} />
        </label>
      </div>

      {/* Error */}
      {error && !loading && (
        <p className="mt-1 text-[11px] text-amber-600">{error}</p>
      )}

      {/* Result cards */}
      {results && results.length > 0 && (
        <div className="mt-1.5 flex max-h-52 flex-col gap-1 overflow-y-auto rounded-lg border border-slate-200 bg-white p-1.5">
          {results.map((hit) => (
            <button key={hit.accession + hit.name} onClick={() => pick(hit)}
              className={`flex w-full items-start gap-2 rounded-md px-2.5 py-1.5 text-left transition-colors hover:bg-teal-50 ${selected === hit.name ? "bg-teal-50 ring-1 ring-teal-400" : ""}`}>
              {/* Name + arm badge */}
              <div className="min-w-0 flex-1">
                <div className="flex items-center gap-1.5 flex-wrap">
                  <span className="font-mono text-[11px] font-semibold text-slate-800">{hit.name}</span>
                  {hit.arm && (
                    <span className={`rounded px-1 py-0.5 text-[9px] font-bold ${ARM_COLOR[hit.arm] || ARM_COLOR[""]}`}>{hit.arm}</span>
                  )}
                  <span className="text-[10px] text-slate-400">{hit.accession}</span>
                </div>
                <div className="mt-0.5 font-mono text-[11px] tracking-wide text-teal-700">{hit.sequence}</div>
                {hit.note && <div className="mt-0.5 text-[10px] text-slate-400 truncate">{hit.note}</div>}
              </div>
              {/* Load indicator */}
              {selected === hit.name
                ? <span className="shrink-0 text-[10px] font-medium text-teal-600 mt-0.5">✓ loaded</span>
                : <span className="shrink-0 text-[10px] text-slate-300 mt-0.5 hover:text-teal-500">load →</span>}
            </button>
          ))}
          <p className="pt-0.5 text-center text-[9px] text-slate-300">
            miRBase 22.1 · Kozomara et al., Nucleic Acids Res. 2019
          </p>
        </div>
      )}
    </div>
  );
}

// Base-level primer alignment: renders the actual input sequence (wrapped) with
// each primer's footprint colored + underlined beneath it. Engine-agnostic — takes
// the same {start,end,strand,label,color} feature shape used by AnnotationMap.
function PrimerAlignment({ seqStr, features = [], snp = null, variants = [], perRow = 60, exons = null, isCdna = false }) {
  if (!seqStr || seqStr.length < 2) return null;
  const seq = seqStr.toUpperCase();
  const N = seq.length;
  const feats = features.filter((f) => f && f.start != null && f.end != null && f.end > f.start);
  const variantPos = new Map(variants.map((v) => [v.pos, v]));

  // Build per-base exon label (null = intron/unknown, exon_number = exonic)
  const exonAtBase = new Array(N).fill(null);
  if (exons) {
    exons.forEach(ex => {
      for (let i = Math.max(0, ex.start); i < Math.min(N, ex.end); i++) {
        exonAtBase[i] = ex.exon_number;
      }
    });
  }

  const byBase = new Array(N).fill(null);
  [...feats].sort((a, b) => (a.end - a.start) - (b.end - b.start)).forEach((f) => {
    for (let i = Math.max(0, f.start); i < Math.min(N, f.end); i++) if (!byBase[i]) byBase[i] = f;
  });
  const rows = [];
  for (let s = 0; s < N; s += perRow) rows.push([s, Math.min(N, s + perRow)]);

  // Build exon label segments for each row (for row-level exon labels)
  const getRowExonSegments = (rs, re) => {
    const segs = [];
    let curExon = exonAtBase[rs], segStart = rs;
    for (let i = rs + 1; i <= re; i++) {
      const e = i < N ? exonAtBase[i] : -1;
      if (e !== curExon) {
        segs.push({ exon: curExon, start: segStart, end: i });
        curExon = e; segStart = i;
      }
    }
    return segs;
  };

  return (
    <div className="mt-3">
      <p className="mb-1 text-[11px] font-semibold text-slate-600">Primer alignment on input sequence{" "}
        {exons && (
          <span className={`ml-1 rounded px-1.5 py-0.5 text-[10px] font-bold ${isCdna ? "bg-teal-100 text-teal-700" : "bg-amber-100 text-amber-700"}`}>
            {isCdna ? "cDNA (spliced, exons contiguous)" : "genomic (introns present)"}
          </span>
        )}
        <span className="ml-1 font-normal text-slate-400">(5′→3′ top strand; colored = primer footprint, underline shows span{exons ? (isCdna ? "; teal = exon, labels mark splice junctions" : "; teal = exon, white = intron") : ""})</span>
      </p>
      <div className="overflow-x-auto rounded border border-slate-200 bg-white p-2">
        {rows.map(([rs, re], ri) => {
          const segs = exons ? getRowExonSegments(rs, re) : [];
          return (
            <div key={ri} className="mb-0.5">
              {/* Exon/intron track — thin colored strip ABOVE the sequence row */}
              {exons && segs.length > 0 && (
                <div className="flex" style={{ height: "4px", marginBottom: "1px" }}>
                  <div className="w-12 shrink-0" /> {/* number gutter spacer */}
                  {segs.map((seg, si) => {
                    const chars = seg.end - seg.start;
                    const pct = (chars / perRow) * 100;
                    return (
                      <div key={si}
                        title={seg.exon != null ? `Exon ${seg.exon}` : "Intron"}
                        style={{ width: `${pct}%` }}
                        className={seg.exon != null ? "bg-teal-400" : "bg-slate-200"} />
                    );
                  })}
                </div>
              )}
              {/* Exon label row — only at exon transitions */}
              {exons && segs.some(s => s.exon != null) && (
                <div className="font-mono text-[9px] leading-3" style={{ whiteSpace: "pre" }}>
                  <span className="mr-2 inline-block w-10 shrink-0" />
                  {segs.map((seg, si) => {
                    const chars = seg.end - seg.start;
                    // Only show label at segment start if it's an exon
                    const label = seg.exon != null
                      ? `E${seg.exon}`.padEnd(chars).slice(0, chars)
                      : " ".repeat(chars);
                    return (
                      <span key={si} className={seg.exon != null ? "text-teal-600 font-bold" : "text-slate-200"}>
                        {label}
                      </span>
                    );
                  })}
                </div>
              )}
              {/* Sequence row — clean white background, primers colored on top */}
              <div className="font-mono text-[11px] leading-5" style={{ whiteSpace: "pre" }}>
                <span className="mr-2 inline-block w-10 shrink-0 select-none text-right text-slate-300">{rs + 1}</span>
                {seq.slice(rs, re).split("").map((b, j) => {
                  const i = rs + j;
                  const f = byBase[i];
                  const isSnp = snp != null && i === snp;
                  const v = variantPos.get(i);
                  const style = {
                    // Primers: colored text + solid underline — clearly visible on white
                    color: isSnp ? "#fff" : (f ? f.color : "#334155"),
                    background: isSnp ? "#e11d48" : (v ? "#fde68a" : "transparent"),
                    borderBottom: f ? `2px solid ${f.color}` : "none",
                    fontWeight: f ? "600" : "normal",
                  };
                  const title = isSnp ? "SNP"
                    : v ? `${v.rsid || "variant"} ${v.ref || ""}>${v.alt || ""}${v.maf ? ` MAF ${v.maf}` : ""}`
                    : exonAtBase[i] != null ? `Exon ${exonAtBase[i]}`
                    : f ? f.label
                    : (exons ? "Intron" : "");
                  return <span key={j} style={style} title={title}>{b}</span>;
                })}
              </div>
            </div>
          );
        })}
      </div>
      <div className="mt-1 flex flex-wrap gap-x-3 gap-y-0.5 text-[10px]">
        {feats.map((f, i) => (
          <span key={i} className="flex items-center gap-1">
            <span className="inline-block h-2 w-3 rounded" style={{ background: f.color }} />
            {f.label} <span className="text-slate-400">({f.start + 1}–{f.end}, {f.strand === "rev" ? "rev ◀" : "fwd ▶"})</span>
          </span>
        ))}
        {exons && <span className="flex items-center gap-1"><span className="inline-block h-1.5 w-4 rounded-sm bg-teal-400" /> exon</span>}
        {exons && <span className="flex items-center gap-1"><span className="inline-block h-1.5 w-4 rounded-sm bg-slate-200" /> intron</span>}
        {variants.length > 0 && <span className="flex items-center gap-1"><span className="inline-block h-2 w-3 rounded bg-yellow-200" /> variant</span>}
        {snp != null && <span className="flex items-center gap-1"><span className="inline-block h-2 w-3 rounded bg-rose-600" /> SNP</span>}
      </div>
    </div>
  );
}

function RiskBadge({ risk }) {  const map = {
    high:     { c: "text-rose-700 bg-rose-50 ring-rose-600/20",     t: "High risk" },
    moderate: { c: "text-amber-700 bg-amber-50 ring-amber-600/20",  t: "Moderate" },
    low:      { c: "text-emerald-700 bg-emerald-50 ring-emerald-600/20", t: "Low risk" },
    expected: { c: "text-sky-700 bg-sky-50 ring-sky-600/20",        t: "Expected ✓" },
  };
  const m = map[risk] || map.low;
  return <span className={`inline-flex items-center rounded-full px-2 py-0.5 text-[11px] font-medium ring-1 ring-inset ${m.c}`}>{m.t}</span>;
}

// Parse gene symbol and locus from an NCBI subject title or accession.
// Returns { gene, locus, organism } — all may be empty string.
function parseHitGene(subjectId, subjectTitle) {
  const id = subjectId || "";
  const t = (subjectTitle || "").toLowerCase();

  // Organism: "Homo sapiens" style
  const orgM = (subjectTitle || "").match(/([A-Z][a-z]+ [a-z]+)\s/);
  const organism = orgM ? orgM[1] : "";

  // Accession prefix → locus type
  const accPfx = id.split("_")[0].toUpperCase();
  const locusType =
    ["XM","NM"].includes(accPfx) ? "mRNA" :
    ["XR","NR"].includes(accPfx) ? "ncRNA" :
    ["XP","NP"].includes(accPfx) ? "protein" :
    accPfx === "NC"              ? "chromosome" :
    accPfx === "CP"              ? "genome" :
    accPfx === "NG"              ? "gene-region" : "";

  // Gene name: extract from parentheses like "(BRCA1)" or "(MTHFR)"
  const geneM = (subjectTitle || "").match(/\(([A-Z][A-Z0-9]{1,15})\)/);
  const gene = geneM ? geneM[1] : "";

  return { gene, locus: locusType, organism };
}

// Group hits: same gene → same group (database redundancy); different gene → separate group.
function groupHits(hits) {
  const groups = {}; // gene -> [hits]
  const noGene = [];
  hits.forEach(h => {
    const { gene } = parseHitGene(h.subject_id, h.subject_title);
    if (gene) {
      if (!groups[gene]) groups[gene] = [];
      groups[gene].push(h);
    } else {
      noGene.push(h);
    }
  });
  return { groups, noGene };
}

function SpecificityPanel({ oligos, backendActive, database, onDatabaseChange, isMirna }) {
  const [results, setResults] = useState({});
  const [openName, setOpenName] = useState(null);

  if (!oligos.length) return null;

  const runCheck = async (name, seq) => {
    setResults((r) => ({ ...r, [name]: "loading" }));
    setOpenName(name);
    const res = await backendSpecificity(seq, database, { isMirna });
    setResults((r) => ({ ...r, [name]: res || { error: "No response from backend." } }));
  };

  return (
    <div className="mt-4 rounded-lg border border-slate-200 bg-white p-4">
      <div className="mb-2 flex items-center justify-between">
        <span className="text-xs font-semibold text-slate-700">Specificity (BLAST)</span>
        {backendActive && (
          <div className="flex items-center gap-1.5">
            <span className="text-[10px] text-slate-400">db:</span>
            <input value={database} onChange={(e) => onDatabaseChange(e.target.value)} spellCheck={false}
              className="w-44 rounded-md border border-slate-200 bg-slate-50 px-2 py-0.5 text-[11px] font-mono outline-none focus:border-teal-500"
              placeholder="ncbi / auto / local-db-name" />
          </div>
        )}
      </div>

      {isMirna && (
        <p className="mb-2 rounded-md bg-sky-50 px-3 py-2 text-[11px] text-sky-700">
          <strong>miRNA mode:</strong> hits to miRNA / ncRNA records are{" "}
          <span className="font-semibold">Expected ✓</span> — target + family members, not off-target.
          Only genuine mRNA / genomic hits are <span className="font-semibold text-rose-700">High risk</span>.
          Family discrimination is handled by the panel below.
        </p>
      )}

      {!backendActive ? (
        <p className="text-xs text-slate-500">Connect a FreePrimers backend to BLAST these primers and check for off-target binding sites.</p>
      ) : (
        <div className="space-y-2">
          {oligos.map((o) => {
            const res = results[o.name];
            const isOpen = openName === o.name;
            return (
              <div key={o.name} className="rounded-md border border-slate-100">
                <div className="flex items-center justify-between gap-2 px-3 py-2">
                  <span className="text-xs font-medium text-slate-700">{o.name}</span>
                  <button onClick={() => runCheck(o.name, o.seq)} disabled={res === "loading"}
                    className="inline-flex items-center gap-1.5 rounded-md bg-slate-900 px-2.5 py-1 text-[11px] font-semibold text-white hover:bg-slate-800 disabled:bg-slate-300">
                    {res === "loading" ? <Loader2 size={11} className="animate-spin" /> : <ShieldCheck size={11} />}
                    {res === "loading" ? "Checking…" : "Check"}
                  </button>
                </div>

                {isOpen && res && res !== "loading" && (
                  <div className="border-t border-slate-100 px-3 py-2 space-y-2">
                    {res.error ? (
                      <p className="text-[11px] text-rose-600">{res.error}</p>
                    ) : (
                      <>
                        {res.warnings?.map((w, i) => (
                          <p key={i} className="text-[11px] text-amber-700">{w}</p>
                        ))}

                        {/* Smart interpretation */}
                        <SpecificityInterpretation hits={res.hits} totalHits={res.total_hits}
                          database={res.database} backend={res.backend_used} isMirna={isMirna} />
                      </>
                    )}
                  </div>
                )}
              </div>
            );
          })}
        </div>
      )}
    </div>
  );
}

// Smart hit interpretation — groups hits by gene, distinguishes database redundancy
// (same gene in many records) from genuine off-targets (different genes).
function SpecificityInterpretation({ hits, totalHits, database, backend, isMirna }) {
  const [expanded, setExpanded] = useState(false);
  if (!hits) return null;

  const { groups, noGene } = groupHits(hits);
  const geneNames = Object.keys(groups);

  // Expected hits (miRNA mode)
  const expectedHits = hits.filter(h => h.risk === "expected");
  const realHits = hits.filter(h => h.risk !== "expected");

  // Among real hits, how many unique genes?
  const { groups: realGroups, noGene: realNoGene } = groupHits(realHits);
  const realGenes = Object.keys(realGroups);
  const offTargetGeneCount = realGenes.length + (realNoGene.length > 0 ? 1 : 0);
  const highRisk = realHits.filter(h => h.risk === "high");

  // Verdict
  let verdict, verdictColor, verdictNote;
  if (realHits.length === 0 && expectedHits.length > 0) {
    verdict = "Target confirmed ✓";
    verdictColor = "text-sky-700";
    verdictNote = `${expectedHits.length} hit${expectedHits.length>1?"s":""} all match miRNA/family records — expected.`;
  } else if (highRisk.length === 0) {
    verdict = "Looks specific ✓";
    verdictColor = "text-emerald-700";
    verdictNote = totalHits === 0
      ? "No hits found in this database — primer appears unique."
      : realGenes.length <= 1
        ? `${totalHits} hit${totalHits>1?"s":""} — all map to ${realGenes[0] ? `one locus (${realGenes[0]})` : "one unidentified locus"}. This is typical database redundancy (same gene in multiple records), not off-target amplification.`
        : `${highRisk.length} high-risk off-target hits.`;
  } else if (offTargetGeneCount <= 1) {
    verdict = "Likely specific ✓";
    verdictColor = "text-emerald-600";
    verdictNote = `${totalHits} hit${totalHits>1?"s":""} — all appear to map to ${realGenes[0] ? `the same locus (${realGenes[0]})` : "one locus"}. Multiple hits to the same gene are normal database redundancy (different assemblies, isoforms, organisms). Not a real off-target.`;
  } else {
    verdict = `Off-target risk — ${offTargetGeneCount} distinct loci`;
    verdictColor = "text-rose-700";
    verdictNote = `${highRisk.length} high-risk hit${highRisk.length>1?"s":""} across ${offTargetGeneCount} different loci. Review these — the primer may co-amplify another gene.`;
  }

  return (
    <div className="space-y-2">
      {/* Summary line */}
      <p className="text-[11px] text-slate-400">
        {totalHits} hit{totalHits===1?"":"s"} against '{database}' via {backend} BLAST.
        {isMirna && expectedHits.length > 0 && (
          <span className="ml-1.5 text-sky-600">{expectedHits.length} expected (miRNA/family) · {highRisk.length} off-target high risk</span>
        )}
      </p>

      {/* Verdict badge */}
      <div className={`rounded-md border px-3 py-2 text-[11px] ${
        verdictColor.includes("emerald") ? "border-emerald-200 bg-emerald-50" :
        verdictColor.includes("sky")     ? "border-sky-200 bg-sky-50" :
                                           "border-rose-200 bg-rose-50"
      }`}>
        <span className={`font-semibold ${verdictColor}`}>{verdict}</span>
        <span className="ml-2 text-slate-600">{verdictNote}</span>
      </div>

      {/* Gene grouping summary */}
      {realGenes.length > 0 && (
        <div className="rounded-md bg-slate-50 px-3 py-2 text-[11px]">
          <p className="mb-1 font-medium text-slate-600">Hits by locus:</p>
          {realGenes.map(gene => (
            <div key={gene} className="flex items-center justify-between gap-2 py-0.5">
              <span className="font-mono font-semibold text-slate-700">{gene}</span>
              <span className="text-slate-400">{realGroups[gene].length} record{realGroups[gene].length>1?"s":""}</span>
              <RiskBadge risk={realGroups[gene][0].risk} />
            </div>
          ))}
          {realNoGene.length > 0 && (
            <div className="flex items-center justify-between gap-2 py-0.5">
              <span className="text-slate-400 italic">Unidentified loci</span>
              <span className="text-slate-400">{realNoGene.length}</span>
              <RiskBadge risk={realNoGene[0].risk} />
            </div>
          )}
        </div>
      )}

      {/* Raw hits (collapsible) */}
      <button onClick={() => setExpanded(e => !e)}
        className="text-[10px] text-slate-400 hover:text-slate-600">
        {expanded ? "▲ Hide" : "▼ Show"} raw hits ({hits.length})
      </button>
      {expanded && (
        <div className="space-y-1">
          {hits.slice(0, 12).map((h, i) => (
            <div key={i} className="flex items-center justify-between gap-2 text-[11px]">
              <span className="min-w-0 truncate font-mono text-slate-500">{h.subject_id}{h.subject_title ? ` · ${h.subject_title.slice(0,60)}` : ""}</span>
              <span className="shrink-0 text-slate-400">{h.pct_identity.toFixed(0)}% id, {h.mismatches} mm</span>
              <RiskBadge risk={h.risk} />
            </div>
          ))}
          {hits.length > 12 && <p className="text-[11px] text-slate-400">+ {hits.length - 12} more</p>}
        </div>
      )}
    </div>
  );
}

// ============================================================
//   VALIDATION / COMPARISON PANEL
//   Paste competing primer sets (Primer1, Primer3, ...) and score
//   them on the SAME metrics as our design, via the backend.
// ============================================================
function ComparisonPanel({ sense, snpIndex, allele1, allele2, ourDesign, backendActive }) {
  const [rows, setRows] = useState(null);
  const [busy, setBusy] = useState(false);
  const [err, setErr] = useState("");
  const [draft, setDraft] = useState({ label: "Primer1", fo: "", ro: "", fi: "", ri: "", fi_detects: allele1, ri_detects: allele2, fi_minus2: "", ri_minus2: "" });

  const run = async () => {
    if (!backendActive) { setErr("Connect a backend to score comparisons (uses real thermodynamics + the discrimination model)."); return; }
    setBusy(true); setErr("");
    try {
      const sets = [];
      if (ourDesign && ourDesign.oligos) {
        const byRole = {};
        ourDesign.oligos.forEach((o) => {
          if (o.label.startsWith("Forward inner")) byRole.fi = cleanDNA(o.seq);
          else if (o.label.startsWith("Reverse inner")) byRole.ri = cleanDNA(o.seq);
          else if (o.label.startsWith("Forward outer")) byRole.fo = cleanDNA(o.seq);
          else if (o.label.startsWith("Reverse outer")) byRole.ro = cleanDNA(o.seq);
        });
        if (byRole.fo && byRole.ro && byRole.fi && byRole.ri) {
          sets.push({ label: "FreePrimers (ours)", fo: byRole.fo, ro: byRole.ro, fi: byRole.fi, ri: byRole.ri,
            fi_detects: ourDesign.arms?.a1 || allele1, ri_detects: ourDesign.arms?.a2 || allele2, fi_minus2: "", ri_minus2: "" });
        }
      }
      // validate the pasted competitor set before sending
      const d = draft;
      if (!d.fo || !d.ro || !d.fi || !d.ri) {
        setErr("Fill in all four primers (FO, RO, FI, RI) before scoring."); setBusy(false); return;
      }
      sets.push({
        label: d.label || "Other tool",
        fo: cleanDNA(d.fo), ro: cleanDNA(d.ro), fi: cleanDNA(d.fi), ri: cleanDNA(d.ri),
        fi_detects: (d.fi_detects || allele1).toUpperCase(),
        ri_detects: (d.ri_detects || allele2).toUpperCase(),
        fi_minus2: d.fi_minus2 || "", ri_minus2: d.ri_minus2 || "",
      });
      const res = await fetchWithTimeout(`${getBackendUrl()}/engines/arms/compare`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sense, snp_index: snpIndex, allele1, allele2, sets }),
      }, 15000);
      if (!res.ok) {
        const t = await res.text();
        setErr(`Backend returned ${res.status}: ${t.slice(0, 200)}`);
        setBusy(false); return;
      }
      const j = await res.json();
      setRows(j.rows || []);
    } catch (e) {
      setErr("Comparison failed: " + (e?.message || "backend unreachable or invalid primers."));
    } finally { setBusy(false); }
  };

  const f = (k, v) => setDraft({ ...draft, [k]: v });
  const cell = "rounded border border-slate-200 bg-slate-50 px-2 py-1 font-mono text-[11px] outline-none focus:border-teal-500";

  return (
    <div className="mt-6 rounded-xl border border-slate-200 bg-white p-5 shadow-sm">
      <div className="flex items-center gap-2 text-slate-700">
        <Sparkles size={15} className="text-teal-600" />
        <h3 className="text-sm font-semibold">Validation panel — compare against another tool</h3>
      </div>
      <p className="mt-1 text-xs text-slate-500">
        Run this SNP through Primer1 or Primer3 yourself, paste their four primers below, and score them on the
        same neutral metrics as our design. Our design is included automatically as the reference row.
      </p>
      <div className="mt-3 grid grid-cols-2 gap-2 sm:grid-cols-4">
        <input className={cell} placeholder="label (e.g. Primer1)" value={draft.label} onChange={(e) => f("label", e.target.value)} />
        <input className={cell} placeholder="FO (fwd outer)" value={draft.fo} onChange={(e) => f("fo", e.target.value)} />
        <input className={cell} placeholder="RO (rev outer)" value={draft.ro} onChange={(e) => f("ro", e.target.value)} />
        <div />
        <input className={cell} placeholder="FI (fwd inner)" value={draft.fi} onChange={(e) => f("fi", e.target.value)} />
        <input className={cell} placeholder="RI (rev inner)" value={draft.ri} onChange={(e) => f("ri", e.target.value)} />
        <input className={cell} placeholder={`FI detects (${allele1}/${allele2})`} value={draft.fi_detects} onChange={(e) => f("fi_detects", e.target.value)} />
        <input className={cell} placeholder={`RI detects (${allele1}/${allele2})`} value={draft.ri_detects} onChange={(e) => f("ri_detects", e.target.value)} />
      </div>
      <button onClick={run} disabled={busy || !backendActive}
        className="mt-3 inline-flex items-center gap-1.5 rounded-md bg-slate-900 px-3 py-1.5 text-xs font-semibold text-white hover:bg-slate-800 disabled:bg-slate-300">
        {busy ? <Loader2 size={12} className="animate-spin" /> : <Sparkles size={12} />} Score comparison
      </button>
      {!backendActive && <p className="mt-2 text-[11px] text-slate-400">Needs a connected backend.</p>}
      {err && <p className="mt-2 text-[11px] text-rose-600">{err}</p>}
      {rows && rows.length > 0 && (
        <div className="mt-4 overflow-x-auto">
          <table className="w-full text-left text-[11px]">
            <thead className="text-slate-400">
              <tr>
                <th className="py-1 pr-3 font-medium">Design</th>
                <th className="py-1 pr-3 font-medium">Tm spread ↓</th>
                <th className="py-1 pr-3 font-medium">Discrim. ↑</th>
                <th className="py-1 pr-3 font-medium">Band sep ↑</th>
                <th className="py-1 pr-3 font-medium">Worst dimer ↑</th>
                <th className="py-1 pr-3 font-medium">Worst hairpin ↑</th>
                <th className="py-1 font-medium">In range</th>
              </tr>
            </thead>
            <tbody className="font-mono">
              {rows.map((r, i) => (
                <tr key={i} className={`border-t border-slate-100 ${r.label.includes("ours") ? "bg-teal-50/50" : ""}`}>
                  <td className="py-1.5 pr-3 font-sans font-medium text-slate-700">{r.label}</td>
                  <td className="py-1.5 pr-3">{r.tm_spread}°C</td>
                  <td className="py-1.5 pr-3">{r.min_disc}</td>
                  <td className="py-1.5 pr-3">{r.band_sep ?? "—"}</td>
                  <td className="py-1.5 pr-3">{r.worst_dimer}</td>
                  <td className="py-1.5 pr-3">{r.worst_hairpin}</td>
                  <td className="py-1.5">{r.in_range ? "✓" : "✗"}</td>
                </tr>
              ))}
            </tbody>
          </table>
          <p className="mt-2 text-[10px] text-slate-400">
            Scored identically for every row via the backend. ↓ lower better, ↑ higher better.
          </p>
        </div>
      )}
    </div>
  );
}

// ============================================================
//   miRNA VALIDATION PANEL (family discrimination, backend)
//   Self-contained: paste a family + pick target + paste competing
//   forward primers; scores each on regime + predicted discrimination.
// ============================================================
function MirnaComparePanel({ backendActive }) {
  const [familyText, setFamilyText] = useState("let-7a UGAGGUAGUAGGUUGUAUAGUU\nlet-7c UGAGGUAGUAGGUUGUAUGGUU\nlet-7f UGAGGUAGUAGAUUGUAUAGUU");
  const [target, setTarget] = useState("let-7a");
  const [candText, setCandText] = useState("Other tool  GCGGCTGAGGTAGTAGGTTGT");
  const [rows, setRows] = useState(null);
  const [busy, setBusy] = useState(false);
  const [err, setErr] = useState("");

  const parseFamily = () => {
    const fam = {};
    familyText.split(/\r?\n/).forEach((ln) => {
      const t = ln.trim();
      if (!t) return;
      const sp = t.split(/\s+/);
      if (sp.length >= 2) fam[sp[0]] = sp[1];
    });
    return fam;
  };

  const run = async () => {
    if (!backendActive) { setErr("Connect a backend to score miRNA discrimination."); return; }
    setBusy(true); setErr("");
    try {
      const family = parseFamily();
      const candidates = candText.split(/\r?\n/).map((ln) => {
        const t = ln.trim(); if (!t) return null;
        const m = t.match(/^(.*?)\s+([ACGTUacgtu]+)$/);
        if (m) return { label: m[1].trim(), forward_primer: m[2] };
        return { label: "candidate", forward_primer: t };
      }).filter(Boolean);
      const res = await fetchWithTimeout(`${getBackendUrl()}/engines/mirna/compare`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ target_name: target, family, candidates }),
      }, 15000);
      const j = await res.json();
      setRows(j.rows || []);
    } catch {
      setErr("Comparison failed (backend unreachable or invalid input).");
    } finally { setBusy(false); }
  };

  const cell = "w-full resize-y rounded border border-slate-200 bg-slate-50 p-2 font-mono text-[11px] outline-none focus:border-teal-500";
  const fam = parseFamily();
  const names = Object.keys(fam);

  return (
    <div className="mt-6 rounded-xl border border-slate-200 bg-white p-5 shadow-sm">
      <div className="flex items-center gap-2 text-slate-700">
        <Sparkles size={15} className="text-teal-600" />
        <h3 className="text-sm font-semibold">miRNA discrimination — validation panel</h3>
      </div>
      <p className="mt-1 text-xs text-slate-500">
        Paste a miRNA family (one "name SEQUENCE" per line), choose the target, and paste competing
        forward primers (one "label SEQUENCE" per line). Each is scored on where its 3′ end lands and
        the predicted per-sibling discrimination (strong forward-terminal vs weak RT-window).
      </p>
      <div className="mt-3 grid gap-3 sm:grid-cols-2">
        <div>
          <label className="text-[11px] font-medium text-slate-500">Family</label>
          <textarea rows={4} className={cell} value={familyText} onChange={(e) => setFamilyText(e.target.value)} />
        </div>
        <div>
          <label className="text-[11px] font-medium text-slate-500">Candidate forward primers</label>
          <textarea rows={4} className={cell} value={candText} onChange={(e) => setCandText(e.target.value)} />
        </div>
      </div>
      <div className="mt-2 flex items-center gap-2">
        <label className="text-[11px] font-medium text-slate-500">Target</label>
        <select value={target} onChange={(e) => setTarget(e.target.value)}
          className="rounded border border-slate-200 bg-slate-50 px-2 py-1 text-[11px]">
          {names.map((n) => <option key={n} value={n}>{n}</option>)}
        </select>
        <button onClick={run} disabled={busy || !backendActive}
          className="ml-auto inline-flex items-center gap-1.5 rounded-md bg-slate-900 px-3 py-1.5 text-xs font-semibold text-white hover:bg-slate-800 disabled:bg-slate-300">
          {busy ? <Loader2 size={12} className="animate-spin" /> : <Sparkles size={12} />} Score
        </button>
      </div>
      {!backendActive && <p className="mt-2 text-[11px] text-slate-400">Needs a connected backend.</p>}
      {err && <p className="mt-2 text-[11px] text-rose-600">{err}</p>}
      {rows && rows.length > 0 && (
        <div className="mt-4 space-y-3">
          {rows.map((r, i) => (
            <div key={i} className="rounded-lg border border-slate-100 p-3">
              <div className="flex items-center justify-between">
                <span className="text-xs font-semibold text-slate-700">{r.label}</span>
                {r.note ? <span className="text-[11px] text-amber-600">{r.note}</span> :
                  <span className="text-[11px] text-slate-500">3′ end at miRNA pos {(r.forward_3p_pos ?? 0) + 1} · min discrimination {r.min_discrimination} log</span>}
              </div>
              {r.siblings && r.siblings.length > 0 && (
                <div className="mt-2 flex flex-wrap gap-2">
                  {r.siblings.map((s, j) => (
                    <span key={j} className={`rounded px-2 py-0.5 text-[10px] font-medium ring-1 ring-inset ${s.regime === "forward-terminal" ? "bg-emerald-50 text-emerald-700 ring-emerald-600/20" : s.regime === "rt-window" ? "bg-amber-50 text-amber-700 ring-amber-600/20" : "bg-slate-100 text-slate-500 ring-slate-200"}`}>
                      vs {s.sibling}: {s.regime === "forward-terminal" ? "STRONG" : s.regime === "rt-window" ? "weak" : "none"} {s.discrimination_log}
                    </span>
                  ))}
                </div>
              )}
            </div>
          ))}
        </div>
      )}
    </div>
  );
}

// ============================================================
//   SNP-AWARE PLACEMENT PANEL (backend)
//   Region + known variants -> ranked safe primer placements.
// ============================================================
function SnpAwarePanel({ backendActive }) {
  const [region, setRegion] = useState("");
  const [variantsText, setVariantsText] = useState("");
  const [offset, setOffset] = useState(0);
  const [primerLen, setPrimerLen] = useState(20);
  const [strand, setStrand] = useState("fwd");
  const [conservative, setConservative] = useState(false);
  const [rows, setRows] = useState(null);
  const [busy, setBusy] = useState(false);
  const [err, setErr] = useState("");

  const parseVariants = () => {
    // accept JSON array, or "pos ref alt [maf] [rsid]" lines.
    // The maf field is optional: if the 4th token isn't a number, treat it as
    // the rsID (so "260 C T rs6265" works as well as "260 C T 0.3 rs6265").
    const t = variantsText.trim();
    if (!t) return [];
    if (t.startsWith("[")) { try { return JSON.parse(t); } catch { return []; } }
    return t.split(/\r?\n/).map((ln) => {
      const s = ln.trim(); if (!s) return null;
      const p = s.split(/[\s,]+/);
      const pos = parseInt(p[0], 10);
      let maf = 0.0, rsid = "";
      if (p[3] !== undefined) {
        const asNum = parseFloat(p[3]);
        if (!isNaN(asNum) && /^[\d.eE+-]+$/.test(p[3])) { maf = asNum; rsid = p[4] || ""; }
        else { rsid = p[3]; }   // 4th token is an rsID, not a maf
      }
      return { pos, ref: (p[1] || "N").toUpperCase(), alt: (p[2] || "N").toUpperCase(), maf, rsid };
    }).filter((v) => v && !isNaN(v.pos));
  };

  const run = async () => {
    if (!backendActive) { setErr("Connect a backend to scan placements."); return; }
    const seq = region.toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, "");
    if (seq.length < primerLen) { setErr("Region is shorter than the primer length."); return; }
    setBusy(true); setErr("");
    try {
      const res = await fetchWithTimeout(`${getBackendUrl()}/engines/snp/scan`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sense: seq, variants: parseVariants(), primer_len: primerLen,
          strand, region_offset: offset, top: 10, conservative }),
      }, 15000);
      if (!res.ok) {
        const t = await res.text();
        setErr(`Backend returned ${res.status}: ${t.slice(0, 200)}`);
        setBusy(false); return;
      }
      const j = await res.json();
      setRows(j.placements || []);
      if ((j.placements || []).length === 0) setErr("No placements returned (check region length and primer length).");
    } catch (e) {
      setErr("Scan failed: " + (e?.message || "backend unreachable or invalid input."));
    } finally { setBusy(false); }
  };

  const onUpload = (e) => {
    const f = e.target.files?.[0]; if (!f) return;
    const r = new FileReader();
    r.onload = () => {
      const text = String(r.result || ""); const lines = text.split(/\r?\n/);
      let s = ""; for (const ln of lines) { if (ln.startsWith(">")) { if (s) break; continue; } s += ln.trim(); }
      setRegion(s.toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, ""));
    };
    r.readAsText(f); e.target.value = "";
  };

  const verdictColor = (v) => v === "clean" ? "bg-emerald-50 text-emerald-700 ring-emerald-600/20"
    : v === "caution" ? "bg-amber-50 text-amber-700 ring-amber-600/20"
    : "bg-rose-50 text-rose-700 ring-rose-600/20";

  return (
    <div className="rounded-xl border border-slate-200 bg-white p-5 shadow-sm">
      <p className="mb-3 text-xs text-slate-500">
        Find the safest primer placements in a region given known population variants. Penalizes primers
        overlapping common variants, weighted by allele frequency and distance from the 3′ end. Pull real
        variants locally with <span className="font-mono">fetch_variants.py</span>, then paste them here.
      </p>
      <div className="grid gap-3 sm:grid-cols-2">
        <div>
          <div className="mb-1 flex items-center justify-between">
            <label className="text-[11px] font-medium text-slate-500">Region sequence</label>
            <div className="flex items-center gap-2">
              {region.length > 0 && (
                <button onClick={() => setRegion("")} title="Clear sequence"
                  className="inline-flex items-center gap-1 text-[11px] text-slate-400 hover:text-rose-600">
                  <X size={11} /> Clear
                </button>
              )}
              <label className="inline-flex cursor-pointer items-center gap-1 text-[11px] text-slate-500 hover:text-slate-700">
                <Upload size={11} /> FASTA
                <input type="file" accept=".fa,.fasta,.txt" className="hidden" onChange={onUpload} />
              </label>
            </div>
          </div>
          <textarea rows={4} value={region} onChange={(e) => setRegion(e.target.value)}
            placeholder="Paste region sequence or upload FASTA"
            className="w-full resize-y rounded border border-slate-200 bg-slate-50 p-2 font-mono text-[11px] outline-none focus:border-teal-500" />
        </div>
        <div>
          <label className="text-[11px] font-medium text-slate-500">Variants (JSON, or "pos ref alt maf rsid" per line)</label>
          <textarea rows={4} value={variantsText} onChange={(e) => setVariantsText(e.target.value)}
            placeholder={'32338500 A G 0.3 rs123\n…or paste vars.json contents'}
            className="w-full resize-y rounded border border-slate-200 bg-slate-50 p-2 font-mono text-[11px] outline-none focus:border-teal-500" />
        </div>
      </div>
      <div className="mt-3 flex flex-wrap items-center gap-3 text-[11px]">
        <label className="flex items-center gap-1">offset
          <input type="number" value={offset} onChange={(e) => setOffset(parseInt(e.target.value || "0", 10))}
            className="w-28 rounded border border-slate-200 bg-slate-50 px-1.5 py-0.5 font-mono" /></label>
        <label className="flex items-center gap-1">primer len
          <input type="number" value={primerLen} onChange={(e) => setPrimerLen(parseInt(e.target.value || "20", 10))}
            className="w-16 rounded border border-slate-200 bg-slate-50 px-1.5 py-0.5 font-mono" /></label>
        <label className="flex items-center gap-1">strand
          <select value={strand} onChange={(e) => setStrand(e.target.value)}
            className="rounded border border-slate-200 bg-slate-50 px-1.5 py-0.5"><option value="fwd">fwd</option><option value="rev">rev</option></select></label>
        <label className="flex items-center gap-1">
          <input type="checkbox" checked={conservative} onChange={(e) => setConservative(e.target.checked)} />
          conservative</label>
        <button onClick={run} disabled={busy || !backendActive}
          className="ml-auto inline-flex items-center gap-1.5 rounded-md bg-slate-900 px-3 py-1.5 text-xs font-semibold text-white hover:bg-slate-800 disabled:bg-slate-300">
          {busy ? <Loader2 size={12} className="animate-spin" /> : <Beaker size={12} />} Scan
        </button>
      </div>
      {!backendActive && <p className="mt-2 text-[11px] text-slate-400">Needs a connected backend. (rsID→region variant pull is the local fetch_variants.py workflow.)</p>}
      {err && <p className="mt-2 text-[11px] text-rose-600">{err}</p>}
      {rows && rows.length > 0 && (
        <div className="mt-4 overflow-x-auto">
          <table className="w-full text-left text-[11px]">
            <thead className="text-slate-400"><tr>
              <th className="py-1 pr-3 font-medium">Position</th><th className="py-1 pr-3 font-medium">Sequence</th>
              <th className="py-1 pr-3 font-medium">Risk</th><th className="py-1 pr-3 font-medium">Overlaps</th>
              <th className="py-1 font-medium">Verdict</th></tr></thead>
            <tbody className="font-mono">
              {rows.map((r, i) => (
                <tr key={i} className="border-t border-slate-100">
                  <td className="py-1.5 pr-3">{r.start}–{r.end}</td>
                  <td className="py-1.5 pr-3">{r.seq}</td>
                  <td className="py-1.5 pr-3">{r.total_risk}</td>
                  <td className="py-1.5 pr-3">{r.n_overlaps}</td>
                  <td className="py-1.5"><span className={`rounded px-2 py-0.5 text-[10px] font-medium ring-1 ring-inset ${verdictColor(r.verdict)}`}>{r.verdict}</span></td>
                </tr>
              ))}
            </tbody>
          </table>
          {region && (
            <PrimerAlignment
              seqStr={region.toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, "")}
              variants={parseVariants().map((v) => ({ ...v, pos: v.pos - offset }))}
              features={rows.slice(0, 3).map((r) => ({
                start: r.start, end: r.end, strand: r.strand,
                label: `${r.verdict} (risk ${r.total_risk})`,
                color: r.verdict === "clean" ? "#059669" : r.verdict === "caution" ? "#d97706" : "#e11d48",
              }))} />
          )}
        </div>
      )}
    </div>
  );
}

// ============================================================
//   FIGURES: melt curves, amplicon schematic, CpG map, conversion view
// ============================================================
function MeltCurveFigure({ curves }) {
  if (!curves) return null;
  const T = curves.temperature, W = 460, H = 180, pad = 32;
  const tlo = T[0], thi = T[T.length - 1];
  const sx = (t) => pad + (t - tlo) / (thi - tlo) * (W - 2 * pad);
  const sy = (v) => H - pad - v * (H - 2 * pad);
  const colors = ["#0d9488", "#d97706", "#7c3aed"];
  const names = Object.keys(curves.melt);
  const path = (ys) => ys.map((v, i) => `${i ? "L" : "M"}${sx(T[i]).toFixed(1)},${sy(v).toFixed(1)}`).join(" ");
  return (
    <div className="mt-3">
      <p className="mb-1 text-[11px] font-semibold text-slate-600">Predicted melt curves <span className="font-normal text-slate-400">(normalised fluorescence vs °C)</span></p>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full rounded border border-slate-200 bg-white">
        <line x1={pad} y1={H - pad} x2={W - pad} y2={H - pad} stroke="#cbd5e1" />
        <line x1={pad} y1={pad} x2={pad} y2={H - pad} stroke="#cbd5e1" />
        {names.map((nm, i) => <path key={nm} d={path(curves.melt[nm])} fill="none" stroke={colors[i % 3]} strokeWidth="1.8" />)}
        <text x={pad} y={H - 8} fontSize="9" fill="#64748b">{tlo}°C</text>
        <text x={W - pad - 18} y={H - 8} fontSize="9" fill="#64748b">{thi}°C</text>
      </svg>
      <div className="mt-1 flex flex-wrap gap-3 text-[10px]">
        {names.map((nm, i) => <span key={nm} className="flex items-center gap-1"><span className="inline-block h-2 w-3 rounded" style={{ background: colors[i % 3] }} /> {nm} (Tm {curves.tm[nm]}°C)</span>)}
      </div>
      <p className="mt-1 text-[10px] italic text-slate-400">{curves.model}</p>
    </div>
  );
}

function DiffCurveFigure({ curves }) {
  if (!curves || !curves.difference) return null;
  const T = curves.temperature, W = 460, H = 150, pad = 32;
  const tlo = T[0], thi = T[T.length - 1];
  const sx = (t) => pad + (t - tlo) / (thi - tlo) * (W - 2 * pad);
  const sy = (v) => H - pad - v * (H - 2 * pad);
  const colors = ["#0d9488", "#d97706", "#7c3aed"];
  const names = Object.keys(curves.difference);
  const path = (ys) => ys.map((v, i) => `${i ? "L" : "M"}${sx(T[i]).toFixed(1)},${sy(v).toFixed(1)}`).join(" ");
  return (
    <div className="mt-3">
      <p className="mb-1 text-[11px] font-semibold text-slate-600">Difference plot <span className="font-normal text-slate-400">(−dF/dT; peak = melt temperature)</span></p>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full rounded border border-slate-200 bg-white">
        <line x1={pad} y1={H - pad} x2={W - pad} y2={H - pad} stroke="#cbd5e1" />
        {names.map((nm, i) => <path key={nm} d={path(curves.difference[nm])} fill="none" stroke={colors[i % 3]} strokeWidth="1.8" />)}
      </svg>
    </div>
  );
}

function AmpliconSchematic({ amplicon, snpIndex, regionLen, primers }) {
  if (!amplicon || amplicon.start == null) return null;
  const W = 460, H = 70, pad = 20;
  const L = regionLen || amplicon.end;
  const sx = (p) => pad + (p / L) * (W - 2 * pad);
  return (
    <div className="mt-3">
      <p className="mb-1 text-[11px] font-semibold text-slate-600">Amplicon schematic</p>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full rounded border border-slate-200 bg-white">
        <line x1={pad} y1={H / 2} x2={W - pad} y2={H / 2} stroke="#e2e8f0" strokeWidth="2" />
        <rect x={sx(amplicon.start)} y={H / 2 - 6} width={sx(amplicon.end) - sx(amplicon.start)} height="12" fill="#99f6e4" opacity="0.6" />
        {primers?.map((p, i) => (
          <rect key={i} x={sx(p.start)} y={H / 2 - 5} width={Math.max(3, sx(p.end) - sx(p.start))} height="10" fill={p.strand === "fwd" ? "#0d9488" : "#d97706"} />
        ))}
        {snpIndex != null && (<>
          <line x1={sx(snpIndex)} y1={H / 2 - 16} x2={sx(snpIndex)} y2={H / 2 + 16} stroke="#dc2626" strokeWidth="1.5" />
          <text x={sx(snpIndex) - 6} y={H / 2 - 19} fontSize="9" fill="#dc2626">SNP</text>
        </>)}
        <text x={pad} y={H - 4} fontSize="8" fill="#94a3b8">0</text>
        <text x={W - pad - 14} y={H - 4} fontSize="8" fill="#94a3b8">{L}</text>
      </svg>
    </div>
  );
}

function CpgMap({ cpgs, regionLen, amplicon, primers }) {
  if (!cpgs || !regionLen) return null;
  const W = 460, H = 64, pad = 20;
  const sx = (p) => pad + (p / regionLen) * (W - 2 * pad);
  return (
    <div className="mt-3">
      <p className="mb-1 text-[11px] font-semibold text-slate-600">CpG map <span className="font-normal text-slate-400">({cpgs.length} CpG sites)</span></p>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full rounded border border-slate-200 bg-white">
        <line x1={pad} y1={H / 2} x2={W - pad} y2={H / 2} stroke="#e2e8f0" strokeWidth="2" />
        {amplicon && amplicon.start != null && <rect x={sx(amplicon.start)} y={H / 2 - 7} width={sx(amplicon.end) - sx(amplicon.start)} height="14" fill="#99f6e4" opacity="0.5" />}
        {primers?.map((p, i) => p.start != null && <rect key={i} x={sx(p.start)} y={H / 2 - 5} width={Math.max(3, sx(p.end) - sx(p.start))} height="10" fill={p.strand === "fwd" ? "#0d9488" : "#d97706"} opacity="0.8" />)}
        {cpgs.map((c, i) => <circle key={i} cx={sx(c)} cy={H / 2} r="2.5" fill="#7c3aed" />)}
        <text x={pad} y={H - 4} fontSize="8" fill="#94a3b8">0</text>
        <text x={W - pad - 18} y={H - 4} fontSize="8" fill="#94a3b8">{regionLen}</text>
      </svg>
      <div className="mt-1 flex flex-wrap gap-3 text-[10px] text-slate-500">
        <span className="flex items-center gap-1"><span className="inline-block h-2 w-2 rounded-full bg-violet-600" /> CpG</span>
        <span className="flex items-center gap-1"><span className="inline-block h-2 w-3 bg-teal-300/60" /> amplicon</span>
        <span className="flex items-center gap-1"><span className="inline-block h-2 w-3 bg-teal-600" /> fwd primer</span>
        <span className="flex items-center gap-1"><span className="inline-block h-2 w-3 bg-amber-600" /> rev primer</span>
      </div>
    </div>
  );
}

function ConversionView({ original, unmeth, meth, cpgs, primers }) {
  if (!original || !unmeth) return null;
  const cpgSet = new Set(cpgs || []);
  const COMPL = { A: "T", T: "A", C: "G", G: "C", N: "N" };
  const rc = (s) => s.split("").reverse().map((b) => COMPL[b] || "N").join("");

  const baseCls = (i, b, compareTo) => {
    const changed = compareTo && compareTo[i] !== b;
    if (cpgSet.has(i)) return "bg-violet-200 text-violet-900";
    if (changed) return "text-rose-600";
    return "text-slate-600";
  };
  const cover = {};
  (primers || []).forEach((p) => { if (p.start != null) for (let i = p.start; i < p.end; i++) cover[i] = p; });

  // per-primer alignment block: primer vs its binding site on the converted template.
  // The methylated/unmethylated template determines what the primer actually pairs with.
  // We align against BOTH converted states so the CpG-discriminating base is visible.
  const PrimerAlignment = ({ p }) => {
    if (p.start == null) return null;
    // the template the primer binds: forward primer matches the top (sense) converted
    // strand; reverse primer matches the bottom strand, i.e. revcomp of the top region.
    const blockU = unmeth.slice(p.start, p.end);
    const blockM = meth.slice(p.start, p.end);
    // primer as it would pair: forward primer == top strand; reverse primer is given
    // 5'->3' on the bottom strand, so its complement reads along the top strand.
    const fwd = p.strand === "fwd";
    // For display we show the primer 5'->3' and the matching template region 5'->3'.
    const primerSeq = p.seq;
    const tmplU = fwd ? blockU : rc(blockU);
    const tmplM = fwd ? blockM : rc(blockM);
    const matchBar = (tmpl) => primerSeq.split("").map((b, i) =>
      (tmpl[i] === b ? "|" : tmpl[i] === undefined ? " " : "·")).join("");
    const seqRow = (seq, tmpl) => (
      <span className="break-all font-mono text-[10px] leading-4">
        {seq.split("").map((b, i) => {
          const m = tmpl && tmpl[i] === b;
          const isCpgInPrimer = fwd ? cpgSet.has(p.start + i) : cpgSet.has(p.end - 1 - i);
          return <span key={i} className={isCpgInPrimer ? "bg-violet-200 text-violet-900" : m ? "text-slate-600" : "text-rose-600"}>{b}</span>;
        })}
      </span>
    );
    return (
      <div className="rounded border border-slate-100 bg-white p-2">
        <p className="mb-0.5 text-[10px] font-semibold text-slate-600">{p.role} <span className="font-normal text-slate-400">({fwd ? "matches sense converted strand" : "binds antisense; shown 5'→3'"}, {p.n_cpg} CpG)</span></p>
        <div className="space-y-0">
          <div className="flex gap-2"><span className="w-20 shrink-0 text-[9px] text-slate-400">primer 5'→3'</span>{seqRow(primerSeq, tmplU)}</div>
          <div className="flex gap-2"><span className="w-20 shrink-0 text-[9px] text-slate-400"> </span><span className="break-all font-mono text-[10px] leading-4 text-slate-300">{matchBar(tmplU)}</span></div>
          <div className="flex gap-2"><span className="w-20 shrink-0 text-[9px] text-slate-400">unmeth tmpl</span><span className="break-all font-mono text-[10px] leading-4 text-slate-500">{tmplU}</span></div>
          <div className="flex gap-2"><span className="w-20 shrink-0 text-[9px] text-slate-400">meth tmpl</span>
            <span className="break-all font-mono text-[10px] leading-4">
              {tmplM.split("").map((b, i) => <span key={i} className={tmplU[i] !== b ? "bg-violet-200 text-violet-900" : "text-slate-500"}>{b}</span>)}
            </span>
          </div>
        </div>
      </div>
    );
  };

  // Block-wrapped alignment: split into fixed-width chunks so the template rows
  // and the primer-footprint rows stay column-aligned even when wrapped.
  const BLOCK = 60;
  const nBlocks = Math.ceil(original.length / BLOCK);
  const monoCell = { display: "inline-block", width: "0.62em", textAlign: "center" };
  const cellSpan = (b, cls, key) => <span key={key} style={monoCell} className={cls}>{b}</span>;

  const block = (bi) => {
    const s = bi * BLOCK, e = Math.min(original.length, s + BLOCK);
    const idxs = [];
    for (let i = s; i < e; i++) idxs.push(i);
    const seqRow = (label, seq, compareTo, labelCls = "text-slate-500") => (
      <div className="whitespace-nowrap">
        <span className={`mr-2 inline-block w-24 text-right text-[10px] font-medium ${labelCls}`}>{label}</span>
        {idxs.map((i) => cellSpan(seq[i], baseCls(i, seq[i], compareTo), i))}
      </div>
    );
    // primer footprint rows: one row per primer in this block, arrows only where it sits
    const primerRows = (primers || []).filter((p) => p.start != null && p.end > s && p.start < e).map((p, pi) => (
      <div key={"p" + pi} className="whitespace-nowrap">
        <span className="mr-2 inline-block w-24 truncate text-right text-[9px] text-slate-400" title={p.role}>{p.role}</span>
        {idxs.map((i) => {
          const inside = i >= p.start && i < p.end;
          if (!inside) return <span key={i} style={monoCell}> </span>;
          const color = p.strand === "fwd" ? "text-teal-700" : "text-amber-700";
          return <span key={i} style={monoCell} className={`${color} font-bold`}>{p.strand === "fwd" ? "▸" : "◂"}</span>;
        })}
      </div>
    ));
    return (
      <div key={bi} className="mb-2 text-[10px] leading-4">
        <div className="mb-0.5 text-[9px] text-slate-400">{s + 1}–{e}</div>
        {seqRow("original", original, null)}
        {seqRow("unmethylated", unmeth, original)}
        {seqRow("methylated", meth, original)}
        {primerRows}
      </div>
    );
  };

  return (
    <div className="mt-3 space-y-3">
      <div>
        <p className="mb-1 text-[11px] font-semibold text-slate-600">Bisulfite conversion + primer alignment <span className="font-normal text-slate-400">(violet = CpG-C; red = converted C→T; teal ▸ fwd, amber ◂ rev — primers aligned under their binding positions)</span></p>
        <div className="max-h-72 overflow-x-auto overflow-y-auto rounded border border-slate-200 bg-slate-50 p-2 font-mono">
          {Array.from({ length: nBlocks }, (_, bi) => block(bi))}
        </div>
      </div>
      {primers && primers.length > 0 && (
        <div>
          <p className="mb-1 text-[11px] font-semibold text-slate-600">Per-primer alignment <span className="font-normal text-slate-400">(primer vs its binding site; violet = CpG discrimination base, red = mismatch to unmeth template)</span></p>
          <div className="space-y-2">
            {primers.map((p, i) => <PrimerAlignment key={i} p={p} />)}
          </div>
        </div>
      )}
    </div>
  );
}

function CpgIslandFigure({ island, original, cpgs }) {
  if (!island || !island.centers || island.centers.length === 0) return null;
  const cpgSet = new Set(cpgs || []);
  const W = 460, H = 130, pad = 32;
  const L = island.region_length;
  const sx = (p) => pad + (p / L) * (W - 2 * pad);
  const xs = island.centers;
  // two tracks: GC% (0-100) and obs/exp (0-2 clamp)
  const gcY = (v) => (H / 2 - 6) - (v / 100) * (H / 2 - pad);
  const oeY = (v) => (H - pad) - (Math.min(v, 2) / 2) * (H / 2 - pad);
  const gcPath = island.gc.map((v, i) => `${i ? "L" : "M"}${sx(xs[i]).toFixed(1)},${gcY(v).toFixed(1)}`).join(" ");
  const oePath = island.oe.map((v, i) => `${i ? "L" : "M"}${sx(xs[i]).toFixed(1)},${oeY(v).toFixed(1)}`).join(" ");
  return (
    <div className="mt-3">
      <p className="mb-1 text-[11px] font-semibold text-slate-600">Predicted CpG island <span className="font-normal text-slate-400">(Gardiner-Garden: GC% &gt; 50 & obs/exp &gt; 0.6)</span></p>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full rounded border border-slate-200 bg-white">
        {island.intervals.map(([s, e], i) => (
          <rect key={i} x={sx(s)} y={pad - 6} width={sx(e) - sx(s)} height={H - pad - (pad - 6)} fill="#a7f3d0" opacity="0.45" />
        ))}
        {/* GC% threshold line (50%) */}
        <line x1={pad} y1={gcY(50)} x2={W - pad} y2={gcY(50)} stroke="#94a3b8" strokeDasharray="3,3" strokeWidth="0.6" />
        {/* obs/exp threshold (0.6) */}
        <line x1={pad} y1={oeY(0.6)} x2={W - pad} y2={oeY(0.6)} stroke="#94a3b8" strokeDasharray="3,3" strokeWidth="0.6" />
        <path d={gcPath} fill="none" stroke="#0d9488" strokeWidth="1.4" />
        <path d={oePath} fill="none" stroke="#7c3aed" strokeWidth="1.4" />
        <text x={pad} y={pad - 8} fontSize="8" fill="#0d9488">GC%</text>
        <text x={pad} y={H / 2 + 6} fontSize="8" fill="#7c3aed">obs/exp CpG</text>
      </svg>
      <p className="mt-1 text-[10px] text-slate-400">
        {island.intervals.length > 0
          ? `Predicted CpG island region(s): ${island.intervals.map(([s, e]) => `${s}–${e}`).join(", ")}.`
          : "No CpG island predicted in this region by the GG&F criteria."}
      </p>
      {original && (
        <div className="mt-2">
          <p className="mb-1 text-[10px] font-medium text-slate-500">Sequence <span className="font-normal text-slate-400">(violet = CpG dinucleotide; green shading = predicted island)</span></p>
          <div className="max-h-32 overflow-y-auto rounded border border-slate-200 bg-slate-50 p-2">
            <span className="break-all font-mono text-[10px] leading-4">
              {original.split("").map((b, i) => {
                const inIsland = island.intervals.some(([s, e]) => i >= s && i < e);
                const isCpg = cpgSet.has(i) || (i > 0 && cpgSet.has(i - 1)); // C or G of a CpG
                let cls = "text-slate-500";
                if (isCpg) cls = "bg-violet-300 font-bold text-violet-900";
                else if (inIsland) cls = "bg-emerald-100 text-emerald-800";
                return <span key={i} className={cls}>{b}</span>;
              })}
            </span>
          </div>
        </div>
      )}
    </div>
  );
}

function DetailsPanel({ title, rows }) {
  const items = rows.filter((r) => r[1] != null && r[1] !== "");
  if (items.length === 0) return null;
  return (
    <div className="mt-3 rounded-lg border border-slate-200 bg-slate-50/60 p-3">
      <p className="mb-1.5 text-[11px] font-semibold text-slate-600">{title}</p>
      <div className="grid grid-cols-2 gap-x-4 gap-y-1 sm:grid-cols-3">
        {items.map(([k, v], i) => (
          <div key={i} className="flex flex-col">
            <span className="text-[10px] uppercase tracking-wide text-slate-400">{k}</span>
            <span className="font-mono text-[12px] font-medium text-slate-700">{v}</span>
          </div>
        ))}
      </div>
    </div>
  );
}

// ── Reusable info cards ───────────────────────────────────────────────────────

// Displayed after a miRNA is selected from the database search
function MirnaInfoCard({ hit }) {
  if (!hit) return null;
  const ARM_COLOR = { "5p": "bg-teal-100 text-teal-700", "3p": "bg-violet-100 text-violet-700" };
  return (
    <div className="mt-2 rounded-md border border-teal-200 bg-teal-50 px-3 py-2 text-[11px]">
      <div className="flex flex-wrap items-center gap-2">
        <span className="font-mono font-semibold text-teal-800">{hit.name}</span>
        {hit.arm && (
          <span className={`rounded px-1.5 py-0.5 text-[10px] font-bold ${ARM_COLOR[hit.arm] || "bg-slate-100 text-slate-600"}`}>{hit.arm}</span>
        )}
        <span className="rounded bg-teal-100 px-1.5 py-0.5 text-teal-700">{hit.accession}</span>
        <span className="font-mono tracking-wide text-teal-700">{hit.sequence}</span>
      </div>
      <div className="mt-1.5 flex flex-wrap gap-x-3 gap-y-0.5 text-teal-700">
        {hit.family && (
          <span><span className="font-medium text-teal-500">Family</span>{" "}{hit.family}</span>
        )}
        <span><span className="font-medium text-teal-500">Length</span>{" "}{hit.sequence?.length} nt</span>
        <span><span className="font-medium text-teal-500">Source</span>{" "}miRBase 22.1</span>
      </div>
      {hit.note && (
        <p className="mt-1 text-[10px] text-teal-600 italic">{hit.note}</p>
      )}
    </div>
  );
}

// Visual exon map — shows exon/intron structure as a colored bar
function ExonMap({ exons, seqLen, strategy, onStrategyChange, selectedExons, onSelectExon }) {
  if (!exons || exons.length === 0) return null;

  const STRATS = [
    { val: "any_exon",    label: "Any exon",
      desc: "Primers anywhere in exonic sequence. Does not prevent gDNA co-amplification.",
      icon: "◇" },
    { val: "junction",    label: "Exon–exon junction",
      desc: "One primer spans a splice junction. Genomic DNA cannot amplify (intron interrupts the primer).",
      icon: "⟨⟩" },
    { val: "neighboring", label: "Neighboring exons",
      desc: "Forward primer on one exon, reverse on the next. gDNA amplicon includes the intron (too large to amplify efficiently).",
      icon: "[ ][ ]" },
  ];

  return (
    <div className="mt-3 rounded-lg border border-slate-200 bg-white p-3">
      {/* Header */}
      <div className="mb-2 flex items-center justify-between">
        <div>
          <span className="text-[12px] font-semibold text-slate-700">Exon / intron map</span>
          <span className="ml-2 text-[11px] text-slate-400">{exons.length} exon{exons.length > 1 ? "s" : ""} · {seqLen?.toLocaleString()} bp region</span>
        </div>
        {selectedExons?.length > 0 && (
          <button onClick={() => onSelectExon && onSelectExon(null)}
            className="text-[10px] text-slate-400 hover:text-rose-500">✕ clear selection</button>
        )}
      </div>

      {/* Visual bar */}
      <div className="relative h-8 w-full overflow-hidden rounded bg-slate-100 border border-slate-200" title="Grey = intron, teal = exon. Click an exon to select it.">
        {exons.map((ex, i) => {
          const left = (ex.start / seqLen) * 100;
          const width = ((ex.end - ex.start) / seqLen) * 100;
          const isSel = selectedExons && selectedExons.includes(ex.exon_number);
          return (
            <div key={i}
              title={`Exon ${ex.exon_number}: pos ${ex.start + 1}–${ex.end} (${ex.length} bp)${onSelectExon ? " — click to select" : ""}`}
              onClick={() => onSelectExon && onSelectExon(ex)}
              style={{ left: `${left}%`, width: `${Math.max(width, 0.6)}%` }}
              className={`absolute top-0 h-full border-r border-white/40 transition-all
                ${onSelectExon ? "cursor-pointer" : ""}
                ${isSel ? "bg-teal-700 z-10" : "bg-teal-500 hover:bg-teal-600"}`}>
              {width > 2.5 && (
                <span className="absolute inset-0 flex items-center justify-center text-[8px] font-bold text-white select-none">
                  {ex.exon_number}
                </span>
              )}
            </div>
          );
        })}
      </div>

      {/* Selected exon info */}
      {selectedExons?.length > 0 && (
        <div className="mt-2 flex flex-wrap gap-2">
          {selectedExons.map(num => {
            const ex = exons.find(e => e.exon_number === num);
            return ex ? (
              <div key={num} className="flex items-center gap-1.5 rounded-md bg-teal-50 border border-teal-300 px-2 py-1 text-[11px]">
                <span className="font-bold text-teal-700">Exon {num}</span>
                <span className="text-teal-600">{ex.length} bp</span>
                <span className="text-teal-500">pos {ex.start + 1}–{ex.end}</span>
                <button onClick={() => onSelectExon && onSelectExon(ex)} className="ml-1 text-slate-400 hover:text-rose-500">✕</button>
              </div>
            ) : null;
          })}
          {strategy === "neighboring" && selectedExons.length === 1 && (
            <div className="flex items-center text-[11px] text-amber-600 italic">← click a second neighboring exon for the reverse primer</div>
          )}
        </div>
      )}

      {/* Compact exon list for large gene (scrollable) */}
      <div className="mt-2 max-h-24 overflow-y-auto">
        <div className="flex flex-wrap gap-1">
          {exons.map((ex, i) => {
            const isSel = selectedExons && selectedExons.includes(ex.exon_number);
            return (
              <button key={i}
                onClick={() => onSelectExon && onSelectExon(ex)}
                title={`Exon ${ex.exon_number}: ${ex.length} bp (pos ${ex.start + 1}–${ex.end})`}
                className={`rounded px-1.5 py-0.5 text-[10px] font-medium border transition-colors
                  ${isSel
                    ? "bg-teal-700 text-white border-teal-700"
                    : "bg-white text-teal-700 border-teal-300 hover:bg-teal-50"}`}>
                E{ex.exon_number}
                <span className="ml-1 text-[9px] opacity-70">{ex.length > 999 ? `${(ex.length/1000).toFixed(1)}k` : ex.length}</span>
              </button>
            );
          })}
        </div>
      </div>

      {/* Strategy selector */}
      {onStrategyChange && (
        <div className="mt-3 border-t border-slate-100 pt-3">
          <p className="mb-2 text-[11px] font-semibold text-slate-700">Primer design strategy
            <span className="ml-1.5 font-normal text-slate-400">— how should primers relate to exon boundaries?</span>
          </p>
          <div className="grid grid-cols-1 gap-2 sm:grid-cols-3">
            {STRATS.map(({ val, label, desc, icon }) => (
              <button key={val} onClick={() => onStrategyChange(val)}
                className={`text-left rounded-lg border p-2.5 text-[11px] transition-all
                  ${strategy === val
                    ? "border-teal-500 bg-teal-50 ring-1 ring-teal-500/40"
                    : "border-slate-200 bg-white hover:border-slate-300"}`}>
                <div className="flex items-center gap-1.5 mb-1">
                  <span className={`font-mono text-[13px] ${strategy === val ? "text-teal-600" : "text-slate-400"}`}>{icon}</span>
                  <span className={`font-semibold ${strategy === val ? "text-teal-800" : "text-slate-700"}`}>{label}</span>
                  {strategy === val && <span className="ml-auto text-teal-500 text-[10px]">✓ active</span>}
                </div>
                <p className="text-[10px] text-slate-500 leading-relaxed">{desc}</p>
                {val === "neighboring" && strategy === "neighboring" && selectedExons?.length === 2 && (
                  <p className="mt-1 text-[10px] text-teal-600 font-medium">E{selectedExons[0]} → E{selectedExons[1]} selected</p>
                )}
              </button>
            ))}
          </div>
          {strategy === "neighboring" && exons.length < 2 && (
            <p className="mt-2 text-[10px] text-rose-600">⚠ Only one exon visible. Fetch a longer region to use neighboring-exon strategy.</p>
          )}
          {strategy === "neighboring" && exons.length >= 2 && selectedExons?.length < 2 && (
            <p className="mt-2 text-[10px] text-amber-600">Select 2 neighboring exons from the bar above — forward primer will target the first, reverse the second.</p>
          )}
          {strategy === "junction" && (
            <p className="mt-2 text-[10px] text-slate-500">Junction-spanning is enforced automatically. If no junction primer can be found, the engine falls back to a warning — switch to Neighboring exons strategy if this happens.</p>
          )}
        </div>
      )}
    </div>
  );
}

// cDNA / transcript fetch for Standard qPCR — fetches the spliced mRNA so exons
// are contiguous and in transcript order. Lets the user pick a transcript.
function CdnaFetchRow({ setSeq, backendActive, onInfo }) {
  const [query, setQuery] = useState("");
  const [status, setStatus] = useState("idle");
  const [note, setNote] = useState("");
  const [transcripts, setTranscripts] = useState([]);
  const [gene, setGene] = useState("");
  const [activeTid, setActiveTid] = useState("");

  const onUpload = (e) => {
    const f = e.target.files?.[0]; if (!f) return;
    const r = new FileReader();
    r.onload = () => {
      const text = String(r.result || ""); let s = "";
      for (const ln of text.split(/\r?\n/)) { if (ln.startsWith(">")) { if (s) break; continue; } s += ln.trim(); }
      setSeq(s.toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, ""));
      setTranscripts([]); if (onInfo) onInfo(null);
    };
    r.readAsText(f); e.target.value = "";
  };

  const listTx = async () => {
    const id = query.trim(); if (!id) return;
    if (!backendActive) { setStatus("failed"); setNote("Connect a backend first (it queries Ensembl)."); return; }
    setStatus("loading"); setNote(""); setTranscripts([]);
    try {
      const res = await fetchWithTimeout(`${getBackendUrl()}/engines/util/list-transcripts`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ query: id }),
      }, 25000);
      const j = await res.json();
      if (j.ok && j.transcripts?.length) {
        setGene(j.gene); setTranscripts(j.transcripts);
        await fetchTx(j.transcripts[0].transcript_id, id);
      } else {
        setStatus("failed"); setNote(j.note || "Gene not found."); if (onInfo) onInfo(null);
      }
    } catch { setStatus("failed"); setNote("Lookup failed (backend unreachable)."); }
  };

  const fetchTx = async (tid, geneQuery) => {
    setStatus("loading"); setNote(""); setActiveTid(tid);
    try {
      const res = await fetchWithTimeout(`${getBackendUrl()}/engines/util/fetch-cdna`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ query: geneQuery || query.trim(), transcript_id: tid }),
      }, 25000);
      const j = await res.json();
      if (j.ok && j.sequence) {
        setSeq(j.sequence);
        setStatus("ok"); setNote("");
        if (onInfo) onInfo({
          type: "cdna", gene: j.gene, assembly: j.assembly, transcriptId: j.transcript_id,
          biotype: j.biotype, seqLen: j.sequence.length, exons: j.exons || null, note: j.note || "",
        });
      } else { setStatus("failed"); setNote(j.note || "cDNA fetch failed."); }
    } catch { setStatus("failed"); setNote("cDNA fetch failed (backend unreachable)."); }
  };

  return (
    <div className="mb-2">
      <div className="flex flex-wrap items-center gap-2">
        <label className="inline-flex cursor-pointer items-center gap-1.5 rounded-md bg-slate-100 px-2.5 py-1 text-[11px] font-medium text-slate-600 hover:bg-slate-200">
          <Upload size={12} /> Upload FASTA
          <input type="file" accept=".fa,.fasta,.txt,.seq" className="hidden" onChange={onUpload} />
        </label>
        <div className="inline-flex items-center gap-1 rounded-md border border-slate-200 bg-slate-50 px-1.5 py-0.5">
          <input value={query} onChange={(e) => setQuery(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && listTx()}
            placeholder="gene symbol (e.g. BRCA2, GAPDH)" spellCheck={false}
            className="w-52 bg-transparent px-1 py-0.5 font-mono text-[11px] outline-none" />
          <button onClick={listTx} disabled={status === "loading" || !backendActive}
            className="inline-flex items-center gap-1 rounded bg-slate-900 px-2 py-1 text-[11px] font-semibold text-white hover:bg-slate-800 disabled:bg-slate-300">
            {status === "loading" ? <Loader2 size={11} className="animate-spin" /> : <Download size={11} />} Fetch cDNA
          </button>
        </div>
      </div>
      <p className="mt-1 text-[10px] text-slate-400">
        Fetches the spliced transcript (exons only, in order) — the correct template for cDNA qPCR. Introns are removed, so junction-spanning primers are designed precisely.
      </p>

      {transcripts.length > 0 && (
        <div className="mt-2 rounded-md border border-slate-200 bg-white p-2">
          <p className="mb-1.5 text-[11px] font-semibold text-slate-600">
            {gene} — {transcripts.length} transcript{transcripts.length > 1 ? "s" : ""}
            <span className="ml-1.5 font-normal text-slate-400">choose one:</span>
          </p>
          <div className="flex flex-col gap-1 max-h-40 overflow-y-auto">
            {transcripts.map((t) => (
              <button key={t.transcript_id}
                onClick={() => fetchTx(t.transcript_id)}
                className={`flex items-center justify-between gap-2 rounded px-2 py-1 text-left text-[11px] transition-colors ${
                  activeTid === t.transcript_id ? "bg-teal-50 ring-1 ring-teal-400" : "hover:bg-slate-50"
                }`}>
                <span className="flex items-center gap-1.5">
                  <span className="font-mono font-medium text-slate-700">{t.transcript_id}</span>
                  {t.is_canonical && <span className="rounded bg-teal-100 px-1 text-[9px] font-bold text-teal-700">canonical</span>}
                  {t.biotype !== "protein_coding" && <span className="rounded bg-slate-100 px-1 text-[9px] text-slate-500">{t.biotype}</span>}
                </span>
                <span className="shrink-0 text-slate-400">{t.length.toLocaleString()} bp · {t.n_exons} exons</span>
              </button>
            ))}
          </div>
        </div>
      )}

      {note && <p className={`mt-1 text-[11px] ${status === "ok" ? "text-teal-700" : "text-amber-600"}`}>{note}</p>}
    </div>
  );
}

// Displayed after a gene/region or rsID is fetched
function SequenceInfoCard({ info, exonStrategy, onExonStrategyChange, selectedExons, onSelectExon }) {
  if (!info) return null;
  // info shape: { type: "region"|"rsid", ...fields }
  if (info.type === "rsid") {
    return (
      <div className="mt-2 rounded-md border border-teal-200 bg-teal-50 px-3 py-2 text-[11px]">
        <div className="flex flex-wrap items-center gap-2 font-semibold text-teal-800">
          <span className="font-mono">{info.rsid}</span>
          <span className="rounded bg-teal-100 px-1.5 py-0.5 font-mono font-bold tracking-wide">
            {info.allele1}/{info.allele2}
          </span>
          {info.assembly && (
            <span className="rounded bg-teal-100 px-1.5 py-0.5 text-teal-700">{info.assembly}</span>
          )}
        </div>
        <div className="mt-1.5 flex flex-wrap gap-x-3 gap-y-0.5 text-teal-700">
          {info.chromosome && (
            <span>
              <span className="font-medium text-teal-500">Chr</span>{" "}
              {info.chromosome}{info.position ? `:${info.position.toLocaleString()}` : ""}
            </span>
          )}
          {info.gene && (
            <span>
              <span className="font-medium text-teal-500">Gene</span>{" "}
              <span className="font-mono font-semibold">{info.gene}</span>
            </span>
          )}
          {info.varClass && (
            <span><span className="font-medium text-teal-500">Type</span>{" "}{info.varClass}</span>
          )}
          {info.consequence && (
            <span><span className="font-medium text-teal-500">Effect</span>{" "}{info.consequence.replace(/_/g, " ")}</span>
          )}
          {info.maf != null && info.maf > 0 && (
            <span><span className="font-medium text-teal-500">MAF</span>{" "}{(info.maf * 100).toFixed(1)}%</span>
          )}
          {info.ancestral && (
            <span><span className="font-medium text-teal-500">Ancestral</span>{" "}{info.ancestral}</span>
          )}
        </div>
        {info.note && <p className="mt-1 text-[10px] text-amber-600">{info.note}</p>}
      </div>
    );
  }
  // region or cdna fetch
  const isCdna = info.type === "cdna";
  return (
    <div className="mt-2 rounded-md border border-teal-200 bg-teal-50 px-3 py-2 text-[11px]">
      <div className="flex flex-wrap items-center gap-2 font-semibold text-teal-800">
        {info.gene && <span className="font-mono">{info.gene}</span>}
        {isCdna && <span className="rounded bg-teal-100 px-1.5 py-0.5 text-[10px] text-teal-700">cDNA / spliced mRNA</span>}
        {info.assembly && (
          <span className="rounded bg-teal-100 px-1.5 py-0.5 text-teal-700">{info.assembly}</span>
        )}
      </div>
      <div className="mt-1.5 flex flex-wrap gap-x-3 gap-y-0.5 text-teal-700">
        {info.region && (
          <span><span className="font-medium text-teal-500">Region</span>{" "}<span className="font-mono">{info.region}</span></span>
        )}
        {info.seqLen && (
          <span><span className="font-medium text-teal-500">Length</span>{" "}{info.seqLen.toLocaleString()} {isCdna ? "nt cDNA" : "bp"}</span>
        )}
        {info.transcriptId && (
          <span><span className="font-medium text-teal-500">Transcript</span>{" "}<span className="font-mono text-[10px]">{info.transcriptId}</span></span>
        )}
        {info.exons && (
          <span><span className="font-medium text-teal-500">Exons</span>{" "}{info.exons.length}</span>
        )}
      </div>
      {info.note && <p className="mt-1 text-[10px] text-amber-600">{info.note}</p>}
      <ExonMap exons={info.exons} seqLen={info.seqLen}
        strategy={exonStrategy} onStrategyChange={onExonStrategyChange}
        selectedExons={selectedExons} onSelectExon={onSelectExon} />
    </div>
  );
}

// shared FASTA-upload + lookup row. lookupMode "rsid" (SNP) or "region" (gene/coords).
function SeqLoadRow({ seq, setSeq, backendActive, onLoaded, onInfo, lookupMode = "rsid" }) {
  const [query, setQuery] = useState("");
  const [status, setStatus] = useState("idle");
  const [note, setNote] = useState("");
  const onUpload = (e) => {
    const f = e.target.files?.[0]; if (!f) return;
    const r = new FileReader();
    r.onload = () => {
      const text = String(r.result || ""); let s = "";
      for (const ln of text.split(/\r?\n/)) { if (ln.startsWith(">")) { if (s) break; continue; } s += ln.trim(); }
      setSeq(s.toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, ""));
      if (onLoaded) onLoaded();
      if (onInfo) onInfo(null);
    };
    r.readAsText(f); e.target.value = "";
  };
  const doFetch = async () => {
    const id = query.trim(); if (!id) return;
    if (!backendActive) { setStatus("failed"); setNote("Connect a backend first."); return; }
    setStatus("loading"); setNote("");
    try {
      if (lookupMode === "region") {
        const res = await fetchWithTimeout(`${getBackendUrl()}/engines/util/fetch-region`, {
          method: "POST", headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ query: id, max_len: 4000 }),
        }, 25000);
        const j = await res.json();
        if (j.ok && j.sequence) {
          setSeq(j.sequence);
          setStatus("ok"); setNote("");
          if (onLoaded) onLoaded();
          if (onInfo) onInfo({ type: "region", gene: j.gene, region: j.region, assembly: j.assembly, seqLen: j.sequence.length, note: j.note || "", exons: j.exons || null, transcriptId: j.transcript_id || "" });
        } else { setStatus("failed"); setNote(j.note || "Lookup failed."); if (onInfo) onInfo(null); }
      } else {
        const res = await fetchWithTimeout(`${getBackendUrl()}/engines/snp/fetch-rsid`, {
          method: "POST", headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ rsid: id, flank: 250 }),
        }, 20000);
        const j = await res.json();
        if (j.ok && j.sense) {
          setSeq(`${j.sense.slice(0, j.snp_index)}[${j.allele1}/${j.allele2}]${j.sense.slice(j.snp_index + 1)}`);
          setStatus("ok"); setNote("");
          if (onLoaded) onLoaded();
          if (onInfo) onInfo({ type: "rsid", rsid: j.rsid, allele1: j.allele1, allele2: j.allele2, assembly: j.assembly, chromosome: j.chromosome, position: j.position, gene: j.gene, consequence: j.consequence, maf: j.maf, ancestral: j.ancestral, varClass: j.var_class, note: j.note || "" });
        } else { setStatus("failed"); setNote(j.note || "Lookup failed."); if (onInfo) onInfo(null); }
      }
    } catch { setStatus("failed"); setNote("Lookup failed (backend unreachable)."); }
  };
  const placeholder = lookupMode === "region" ? "gene or chr:start-end" : "rsID";
  return (
    <div className="mb-2">
      <div className="flex flex-wrap items-center gap-2">
        <label className="inline-flex cursor-pointer items-center gap-1.5 rounded-md bg-slate-100 px-2.5 py-1 text-[11px] font-medium text-slate-600 hover:bg-slate-200">
          <Upload size={12} /> Upload FASTA
          <input type="file" accept=".fa,.fasta,.txt,.seq" className="hidden" onChange={onUpload} />
        </label>
        <div className="inline-flex items-center gap-1 rounded-md border border-slate-200 bg-slate-50 px-1.5 py-0.5">
          <input value={query} onChange={(e) => setQuery(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && doFetch()}
            placeholder={placeholder} spellCheck={false}
            className={`bg-transparent px-1 py-0.5 font-mono text-[11px] outline-none ${lookupMode === "region" ? "w-44" : "w-28"}`} />
          <button onClick={doFetch} disabled={status === "loading" || !backendActive}
            className="inline-flex items-center gap-1 rounded bg-slate-900 px-2 py-1 text-[11px] font-semibold text-white hover:bg-slate-800 disabled:bg-slate-300">
            {status === "loading" ? <Loader2 size={11} className="animate-spin" /> : <Download size={11} />} Fetch
          </button>
        </div>
      </div>
      {lookupMode === "region" && <p className="mt-1 text-[10px] text-slate-400">Enter a gene symbol (e.g. BRCA2, MGMT), an Ensembl ID, or coordinates like 13:32315000-32316000.</p>}
      {note && <p className={`mt-1 text-[11px] ${status === "ok" ? "text-teal-700" : "text-amber-600"}`}>{note}</p>}
    </div>
  );
}

// ============================================================
//   HRM PANEL
// ============================================================
function HrmPanel({ backendActive, onSaveDesign, justSaved }) {
  const [seq, setSeq] = useState("");
  const [snpIdx, setSnpIdx] = useState(null);
  const [a1, setA1] = useState("C");
  const [a2, setA2] = useState("T");
  const [mode, setMode] = useState("standard");
  const [res, setRes] = useState(null);
  const [busy, setBusy] = useState(false);
  const [err, setErr] = useState("");
  const [specDb, setSpecDb] = useState("ncbi");
  const [hrmInfo, setHrmInfo] = useState(null);

  const clean = seq.toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, "");
  const bracket = seq.match(/\[([ACGT])\/([ACGT])\]/i);

  const run = async () => {
    if (!backendActive) { setErr("Connect a backend to design HRM assays."); return; }
    let s = clean, idx = snpIdx, al1 = a1, al2 = a2;
    if (bracket) {
      const pre = seq.slice(0, seq.indexOf("[")).toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, "");
      idx = pre.length; al1 = bracket[1].toUpperCase(); al2 = bracket[2].toUpperCase();
      s = pre + al1 + seq.slice(seq.indexOf("]") + 1).toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, "");
    }
    if (idx == null) { setErr("Click a base below (or paste a [C/T] bracket) to mark the SNP."); return; }
    setBusy(true); setErr("");
    try {
      const r = await fetchWithTimeout(`${getBackendUrl()}/engines/hrm/design`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sense: s, snp_index: idx, allele1: al1, allele2: al2, mode }),
      }, 15000);
      if (!r.ok) { setErr(`Backend returned ${r.status}: ${(await r.text()).slice(0, 200)}`); setBusy(false); return; }
      const j = await r.json(); j._snpIdx = idx; j._regionLen = s.length; j._seq = s; setRes(j);
    } catch (e) { setErr("HRM design failed: " + (e?.message || "")); }
    finally { setBusy(false); }
  };

  const cls = res?.snp_class;
  const clsColor = cls >= 3 ? "text-amber-700" : "text-emerald-700";
  return (
    <div className="rounded-xl border border-slate-200 bg-white p-5 shadow-sm">
      <p className="mb-3 text-xs text-slate-500">
        HRM genotyping by melt curve. <b>Standard</b>: one flanking pair, genotype by amplicon melt
        (best for class 1–2 SNPs). <b>Allele-specific</b>: two allele-specific primers (3′ base on the SNP) + a common primer in one tube; a 5′ GC tag shifts the two allele products apart in Tm so the melt separates all three genotypes (for class 3–4 SNPs).
      </p>
      <SeqLoadRow seq={seq} setSeq={(v) => { setSeq(v); setSnpIdx(null); }} backendActive={backendActive}
        lookupMode="rsid" onInfo={(info) => setHrmInfo(info)} />
      {hrmInfo && <SequenceInfoCard info={hrmInfo} />}
      <div className="mb-2 flex gap-2">
        {[["standard", "Standard HRM"], ["as", "Allele-specific HRM"]].map(([id, lab]) => (
          <button key={id} onClick={() => setMode(id)}
            className={`rounded-md px-2.5 py-1 text-[11px] font-medium ${mode === id ? "bg-slate-900 text-white" : "bg-slate-100 text-slate-600 hover:bg-slate-200"}`}>{lab}</button>
        ))}
      </div>
      <textarea rows={3} value={seq} onChange={(e) => { setSeq(e.target.value); setSnpIdx(null); }}
        placeholder="Paste region with the SNP, or paste [C/T] bracket form, or upload FASTA / fetch rsID above"
        className="w-full resize-y rounded border border-slate-200 bg-slate-50 p-2 font-mono text-[11px] outline-none focus:border-teal-500" />
      {!bracket && clean.length >= 40 && (
        <div className="mt-2">
          <p className="mb-1 text-[11px] text-slate-500">{snpIdx == null ? "Click the SNP base:" : `SNP at position ${snpIdx + 1} (${clean[snpIdx]}).`}</p>
          <div className="max-h-24 overflow-y-auto rounded border border-slate-200 bg-slate-50 p-2 font-mono text-[11px] leading-5">
            {clean.split("").map((b, i) => (
              <span key={i} onClick={() => { setSnpIdx(i); setA1(b); setA2(b === "C" ? "T" : b === "A" ? "G" : "A"); }}
                className={`cursor-pointer rounded-sm px-[0.5px] ${snpIdx === i ? "bg-teal-500 text-white" : "hover:bg-teal-200"}`}>{b}</span>
            ))}
          </div>
        </div>
      )}
      <div className="mt-2 flex flex-wrap items-center gap-2 text-[11px]">
        {!bracket && snpIdx != null && (<>
          <span className="text-slate-500">alleles</span>
          <input value={a1} onChange={(e) => setA1(e.target.value.toUpperCase().slice(0, 1))} className="w-10 rounded border border-slate-200 bg-slate-50 px-1.5 py-0.5 text-center font-mono" />
          <span>/</span>
          <input value={a2} onChange={(e) => setA2(e.target.value.toUpperCase().slice(0, 1))} className="w-10 rounded border border-slate-200 bg-slate-50 px-1.5 py-0.5 text-center font-mono" />
        </>)}
        <button onClick={run} disabled={busy || !backendActive}
          className="ml-auto inline-flex items-center gap-1.5 rounded-md bg-slate-900 px-3 py-1.5 text-xs font-semibold text-white hover:bg-slate-800 disabled:bg-slate-300">
          {busy ? <Loader2 size={12} className="animate-spin" /> : <FlaskConical size={12} />} Design HRM
        </button>
      </div>
      {!backendActive && <p className="mt-2 text-[11px] text-slate-400">Needs a connected backend.</p>}
      {err && <p className="mt-2 text-[11px] text-rose-600">{err}</p>}
      {res && (
        <div className="mt-4 space-y-2">
          {cls > 0 && <p className={`text-xs font-medium ${clsColor}`}>SNP HRM class {cls} {cls >= 3 ? "— hard for amplicon melt" : "— well-suited"}</p>}
          {res.primers?.map((p, i) => (
            <div key={i} className="rounded border border-slate-100 p-2 font-mono text-[11px]">
              <span className="font-sans font-medium text-slate-700">{p.role}: </span>{p.seq}
              <span className="text-slate-400"> (Tm {p.tm}, GC {p.gc}%)</span>
              {p.note && <div className="font-sans text-[10px] text-slate-400">{p.note}</div>}
            </div>
          ))}
          <DetailsPanel title="Assay details" rows={[
            ["Mode", res.mode === "as-hrm" ? "Allele-specific HRM" : "Standard HRM"],
            ["SNP HRM class", res.snp_class],
            ["Amplicon size", (() => {
              const core = res.amplicon?.size;
              if (!core) return null;
              if (res.mode === "as-hrm") {
                const p1len = res.amplicon?.product_a1_len;
                const p2len = res.amplicon?.product_a2_len;
                if (p1len && p2len) {
                  return `${core} bp core (templated) + ${p1len - core} / ${p2len - core} bp 5′ tag → products ${p1len} / ${p2len} bp`;
                }
              }
              return `${core} bp`;
            })()],
            ["Amplicon span", res.amplicon?.start != null ? `${res.amplicon.start}–${res.amplicon.end}` : null],
            ["SNP in amplicon", res.amplicon?.snp_offset_in_amplicon != null ? `+${res.amplicon.snp_offset_in_amplicon}` : null],
            ["Primers", res.primers?.length],
            ["Predicted Tm range", res.melt_curves ? `${Math.min(...Object.values(res.melt_curves.tm))}–${Math.max(...Object.values(res.melt_curves.tm))} °C` : null],
            ["Allele product Tms", res.mode === "as-hrm" && res.amplicon?.product_a1_tm != null ? `${res.amplicon.product_a1_tm} / ${res.amplicon.product_a2_tm} °C (gap ${res.amplicon.tm_gap} °C)` : null],
            ["5′ GC tags (a1 / a2)", res.mode === "as-hrm" && res.amplicon?.tag_a2 ? `${res.amplicon.tag_a1} / ${res.amplicon.tag_a2}` : null],
          ]} />
          {res.mode === "as-hrm" && res.amplicon?.tag_a2 && (
            <div className="mt-2 rounded-md border border-sky-200 bg-sky-50 px-3 py-2 text-[11px] text-sky-800">
              <p className="font-semibold mb-1">How the 5′ GC tag works — why primers bind specifically despite the mismatches</p>
              <p className="text-sky-700">The tag (<span className="font-mono">{res.amplicon.tag_a1}</span> on allele 1, <span className="font-mono">{res.amplicon.tag_a2}</span> on allele 2) is a <strong>non-templated 5′ tail</strong> — the same principle as M13, T7 promoter, and restriction-site tails used routinely in molecular biology.</p>
              <ul className="mt-1 ml-3 list-disc space-y-0.5 text-sky-700">
                <li><strong>Cycle 1:</strong> Only the 3′ templated portion anneals to the genomic DNA. The 5′ tag dangles freely. The 3′ end determines specificity — including the SNP allele base. Extension produces a new strand that now includes the tag.</li>
                <li><strong>Cycle 2+:</strong> The new strand is the template. The full primer (tag + body) matches it perfectly. From this point the tag is part of every PCR product.</li>
                <li><strong>Result:</strong> The two allele-specific products differ in GC content because of the tag → different melt temperatures → HRM separates all three genotypes.</li>
              </ul>
              <p className="mt-1 text-sky-600 text-[10px]">BLAST will show the tag bases as mismatches to genomic DNA — this is expected and not a design flaw. The amplicon core ({res.amplicon.size} bp) is fully templated; the products are {res.amplicon.product_a1_len} / {res.amplicon.product_a2_len} bp including the tags.</p>
            </div>
          )}
          <AmpliconSchematic amplicon={res.amplicon} snpIndex={res._snpIdx} regionLen={res._regionLen} primers={res.primers} />
          {res._seq && res.primers && (
            <PrimerAlignment seqStr={res._seq} snp={res._snpIdx}
              features={res.primers.map((p, i) => ({
                start: p.start, end: p.end, strand: p.strand,
                label: p.role || `primer ${i + 1}`,
                color: p.strand === "rev" ? "#d97706" : (i === 1 ? "#7c3aed" : "#0d9488"),
              }))} />
          )}
          <MeltCurveFigure curves={res.melt_curves} />
          <DiffCurveFigure curves={res.melt_curves} />
          {res.notes?.map((n, i) => <p key={i} className="text-[11px] text-slate-500">{n}</p>)}
          {res.warnings?.map((w, i) => <p key={i} className="text-[11px] text-amber-600">⚠ {w}</p>)}
          <SpecificityPanel
            oligos={getExportOligos(res)}
            backendActive={backendActive}
            database={specDb}
            onDatabaseChange={setSpecDb}
            isMirna={false}
          />
          {onSaveDesign && (
            <button onClick={() => onSaveDesign(res, "hrm", getExportOligos(res))}
              className={`mt-3 inline-flex items-center gap-1.5 rounded-md px-2.5 py-1.5 text-xs font-semibold transition-colors ${justSaved ? "bg-teal-600 text-white" : "bg-slate-900 text-white hover:bg-slate-800"}`}>
              {justSaved ? <><Check size={12} /> Saved</> : <><Save size={12} /> Save design</>}
            </button>
          )}
        </div>
      )}
    </div>
  );
}
function MethylationPanel({ backendActive, onSaveDesign, justSaved }) {
  const [seq, setSeq] = useState("");
  const [mode, setMode] = useState("msp");
  const [res, setRes] = useState(null);
  const [busy, setBusy] = useState(false);
  const [err, setErr] = useState("");
  const [specDb, setSpecDb] = useState("ncbi");
  const [methylInfo, setMethylInfo] = useState(null);

  const run = async () => {
    if (!backendActive) { setErr("Connect a backend to design methylation assays."); return; }
    const s = seq.toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, "");
    if (s.length < 40) { setErr("Paste a longer region (≥40 bp with CpG sites)."); return; }
    setBusy(true); setErr("");
    try {
      const r = await fetchWithTimeout(`${getBackendUrl()}/engines/methylation/design`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sense: s, mode }),
      }, 15000);
      if (!r.ok) { setErr(`Backend returned ${r.status}: ${(await r.text()).slice(0, 200)}`); setBusy(false); return; }
      const j = await r.json(); j._original = s; setRes(j);
    } catch (e) { setErr("Methylation design failed: " + (e?.message || "")); }
    finally { setBusy(false); }
  };

  const allPrimers = res ? Object.values(res.sets || {}).flat() : [];
  return (
    <div className="rounded-xl border border-slate-200 bg-white p-5 shadow-sm">
      <p className="mb-3 text-xs text-slate-500">
        Bisulfite methylation assays. <b>MSP</b>: methylated + unmethylated primer sets covering CpGs.
        <b> BSP</b>: one CpG-free pair (read by sequencing). <b>MS-HRM</b>: CpG-free primers around a CpG-rich
        amplicon, methylation read by melt shift.
      </p>
      <SeqLoadRow seq={seq} setSeq={setSeq} backendActive={backendActive} lookupMode="region"
        onInfo={(info) => setMethylInfo(info)} />
      {methylInfo && <SequenceInfoCard info={methylInfo} />}
      <div className="mb-2 flex gap-2">
        {[["msp", "MSP"], ["bsp", "BSP"], ["mshrm", "MS-HRM"]].map(([id, lab]) => (
          <button key={id} onClick={() => setMode(id)}
            className={`rounded-md px-2.5 py-1 text-[11px] font-medium ${mode === id ? "bg-slate-900 text-white" : "bg-slate-100 text-slate-600 hover:bg-slate-200"}`}>{lab}</button>
        ))}
      </div>
      <textarea rows={3} value={seq} onChange={(e) => setSeq(e.target.value)}
        placeholder="Paste the genomic region (with CpG sites) — NOT pre-converted; the tool does bisulfite conversion. Or upload FASTA / fetch rsID above."
        className="w-full resize-y rounded border border-slate-200 bg-slate-50 p-2 font-mono text-[11px] outline-none focus:border-teal-500" />
      <button onClick={run} disabled={busy || !backendActive}
        className="mt-2 inline-flex items-center gap-1.5 rounded-md bg-slate-900 px-3 py-1.5 text-xs font-semibold text-white hover:bg-slate-800 disabled:bg-slate-300">
        {busy ? <Loader2 size={12} className="animate-spin" /> : <Beaker size={12} />} Design {mode.toUpperCase()}
      </button>
      {!backendActive && <p className="mt-2 text-[11px] text-slate-400">Needs a connected backend.</p>}
      {err && <p className="mt-2 text-[11px] text-rose-600">{err}</p>}
      {res && (
        <div className="mt-4 space-y-3">
          <p className="text-[11px] text-slate-500">{res.assay} · {res.n_cpg_region} CpG sites in region</p>
          {Object.entries(res.sets || {}).map(([setName, primers]) => (
            <div key={setName}>
              {setName !== "primary" && <p className="text-xs font-semibold capitalize text-slate-700">{setName} set</p>}
              {primers.map((p, i) => (
                <div key={i} className="rounded border border-slate-100 p-2 font-mono text-[11px]">
                  <span className="font-sans font-medium text-slate-700">{p.role}: </span>{p.seq}
                  <span className="text-slate-400"> (Tm {p.tm}, GC {p.gc}%, {p.n_cpg} CpG)</span>
                </div>
              ))}
            </div>
          ))}
          <CpgMap cpgs={res.cpg_positions} regionLen={res.region_length} amplicon={res.amplicon} primers={allPrimers} />
          <CpgIslandFigure island={res.cpg_island} original={res._original} cpgs={res.cpg_positions} />
          {res._original && allPrimers.length > 0 && (
            <PrimerAlignment seqStr={res._original}
              variants={(res.cpg_positions || []).map((p) => ({ pos: p, ref: "CpG" }))}
              features={allPrimers.map((p, i) => ({
                start: p.start, end: p.end, strand: p.strand,
                label: p.role || `primer ${i + 1}`,
                color: p.strand === "rev" ? "#d97706" : "#0d9488",
              }))} />
          )}
          {res.melt_curves && (<>
            <MeltCurveFigure curves={res.melt_curves} />
            <DiffCurveFigure curves={res.melt_curves} />
          </>)}
          <DetailsPanel title="Assay details" rows={[
            ["Assay", res.assay],
            ["Amplicon size", res.amplicon?.size ? `${res.amplicon.size} bp` : null],
            ["Amplicon span", res.amplicon?.start != null ? `${res.amplicon.start}–${res.amplicon.end}` : null],
            ["CpGs in region", res.n_cpg_region],
            ["CpGs in amplicon", res.amplicon?.cpgs_in_amplicon],
            ["Primers", allPrimers.length],
            ["Est. Tm shift (meth−unmeth)", res.amplicon?.est_tm_shift_meth_vs_unmeth != null ? `${res.amplicon.est_tm_shift_meth_vs_unmeth} °C` : null],
          ]} />
          <ConversionView original={res._original} unmeth={res.converted_unmethylated} meth={res.converted_methylated} cpgs={res.cpg_positions} primers={allPrimers} />
          {res.amplicon && <p className="text-[11px] text-slate-500">Amplicon: {JSON.stringify(res.amplicon)}</p>}
          {res.notes?.map((n, i) => <p key={i} className="text-[11px] text-slate-500">{n}</p>)}
          {res.warnings?.map((w, i) => <p key={i} className="text-[11px] text-amber-600">⚠ {w}</p>)}
          <SpecificityPanel
            oligos={getExportOligos(res)}
            backendActive={backendActive}
            database={specDb}
            onDatabaseChange={setSpecDb}
            isMirna={false}
          />
          {onSaveDesign && (
            <button onClick={() => onSaveDesign(res, "methyl", getExportOligos(res))}
              className={`mt-3 inline-flex items-center gap-1.5 rounded-md px-2.5 py-1.5 text-xs font-semibold transition-colors ${justSaved ? "bg-teal-600 text-white" : "bg-slate-900 text-white hover:bg-slate-800"}`}>
              {justSaved ? <><Check size={12} /> Saved</> : <><Save size={12} /> Save design</>}
            </button>
          )}
        </div>
      )}
    </div>
  );
}


export default function FreePrimers() {
  const [tab, setTab] = useState("mirna");
  const [method, setMethod] = useState("stemloop");
  const [seq, setSeq] = useState("");
  const [adv, setAdv] = useState(false);
  const [result, setResult] = useState(null);
  const [selectedPair, setSelectedPair] = useState(0);
  const [hoverKey, setHoverKey] = useState(null);
  const [activeEnzymes, setActiveEnzymes] = useState(["EcoRI", "BamHI", "HindIII"]);
  const [p, setP] = useState({ na: 50, mg: 3, dntp: 0.8, primer: 250, target: 60, amplicon: [80, 150] });
  const [backendUrlInput, setBackendUrlInput] = useState("");
  const [backendStatus, setBackendStatus] = useState("disconnected");
  const [showBackendPanel, setShowBackendPanel] = useState(false);
  const [specDatabase, setSpecDatabase] = useState("ncbi");
  const [rsidInput, setRsidInput] = useState("");
  const [rsidStatus, setRsidStatus] = useState("idle");
  const [rsidNote, setRsidNote] = useState("");
  const [rsidDetails, setRsidDetails] = useState(null);
  // Info cards — shown below the fetch row after any successful fetch/select
  const [seqInfo, setSeqInfo] = useState(null);   // SequenceInfoCard data (region or rsid)
  const [mirnaHit, setMirnaHit] = useState(null); // MirnaInfoCard data (selected miRNA record)
  const { designs: savedDesigns, save: saveDesign, remove: removeDesign, clear: clearDesigns } = useSavedDesigns();
  const [workspaceOpen, setWorkspaceOpen] = useState(false);
  const [justSaved, setJustSaved] = useState(false);
  const [exonStrategy, setExonStrategy] = useState("junction");
  const [selectedExons, setSelectedExons] = useState([]);

  // Fetch flanking sequence for an rsID via the backend (Ensembl). Falls back
  // gracefully when no backend is connected or the lookup fails.
  const fetchRsid = async () => {
    const id = rsidInput.trim();
    if (!id) return;
    const url = getBackendUrl();
    if (!url) { setRsidStatus("failed"); setRsidNote("Connect a backend first (rsID lookup queries Ensembl through it)."); return; }
    setRsidStatus("loading"); setRsidNote("");
    try {
      const res = await fetchWithTimeout(`${url}/engines/snp/fetch-rsid`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ rsid: id, flank: 250 }),
      }, 20000);
      const j = await res.json();
      if (j.ok && j.sense) {
        // Insert the variant in bracket form so designTetraPrimer can parse it.
        const i = j.snp_index;
        const bracket = `${j.sense.slice(0, i)}[${j.allele1}/${j.allele2}]${j.sense.slice(i + 1)}`;
        setSeq(bracket);
        setRsidStatus("ok");
        setRsidNote(j.note || "");
        setRsidDetails({
          rsid: j.rsid, allele1: j.allele1, allele2: j.allele2,
          assembly: j.assembly || "", chromosome: j.chromosome || "",
          position: j.position || null, gene: j.gene || "",
          consequence: j.consequence || "", maf: j.maf ?? null,
          varClass: j.var_class || "",
        });
      } else {
        setRsidStatus("failed");
        setRsidDetails(null);
        setRsidNote(j.error || j.note || "Lookup failed. Upload a FASTA or paste the sequence instead.");
      }
    } catch {
      setRsidStatus("failed");
      setRsidDetails(null);
      setRsidNote("Lookup failed (backend unreachable or timed out). Upload a FASTA or paste the sequence instead.");
    }
  };

  const connectBackend = async () => {
    const url = backendUrlInput.trim();
    if (!url) { setBackendUrl(""); setBackendStatus("disconnected"); return; }
    setBackendStatus("checking");
    const ok = await checkBackendHealth(url);
    if (ok) { setBackendUrl(url); setBackendStatus("connected"); }
    else { setBackendUrl(""); setBackendStatus("failed"); }
  };

  const disconnectBackend = () => {
    setBackendUrl(""); setBackendUrlInput(""); setBackendStatus("disconnected");
  };

  const backendActive = backendStatus === "connected";

  const placeholder = tab === "mirna"
    ? "Paste a mature miRNA, 5′→3′  (e.g. UAGCUUAUCAGACUGAUGUUGA). U or T both fine."
    : tab === "snp"
    ? "Paste the region around your SNP with the variant in brackets, e.g. …GCTACG[A/G]TTACGG… (≥150 nt each side recommended)."
    : "Paste your target cDNA / mRNA sequence, 5′→3′.";

  const run = () => {
    setSelectedPair(0);
    setHoverKey(null);
    try {
      if (tab === "standard") setResult(designStandard(seq, p, {
        exons: seqInfo?.exons || null,
        strategy: exonStrategy,
        selectedExons,
      }));
      else if (tab === "snp") setResult(designTetraPrimer(seq, p));
      else setResult(method === "stemloop" ? designStemLoop(seq, p) : designPolyA(seq, p));
    } catch {
      setResult({ error: "Something went wrong while designing. Check the sequence and try again." });
    }
  };

  const reset = (next) => { setTab(next); setResult(null); setSelectedPair(0); setHoverKey(null); setSeqInfo(null); setMirnaHit(null); };

  // Save the current design to the workspace
  const saveCurrentDesign = () => {
    if (!result || result.error) return;
    const oligos = getExportOligos(result, selectedPair);
    if (!oligos.length) return;
    // Build a human name + summary
    const baseName = seqInfo?.gene || mirnaHit?.name ||
      (tab === "snp" && seqInfo?.rsid) || `${TAB_LABELS[tab]} design`;
    const name = `${baseName}`;
    let summary = "";
    if (result.pairs) summary = `${result.pairs.length} pairs, ${result.pairs[selectedPair]?.amplicon || result.pairs[0]?.amplicon} bp`;
    else if (tab === "mirna") summary = mirnaHit?.accession || "stem-loop assay";
    saveDesign({
      name, tab, summary, oligoCount: oligos.length, oligos,
      // store enough to reload the input
      restore: { seq, selectedPair, seqInfo, mirnaHit, exonStrategy, selectedExons, method },
    });
    setJustSaved(true);
    setTimeout(() => setJustSaved(false), 1600);
  };

  // Save handler for sub-panels (HRM, methylation) which hold their own result state
  const saveExternalDesign = (res, panelTab, oligos) => {
    if (!res || !oligos?.length) return;
    const baseName = res._geneName || res.amplicon?.gene ||
      `${TAB_LABELS[panelTab]} design`;
    const summary = res.amplicon?.size ? `${res.amplicon.size} bp amplicon` : `${oligos.length} oligos`;
    saveDesign({ name: baseName, tab: panelTab, summary, oligoCount: oligos.length, oligos, restore: {} });
    setJustSaved(true);
    setTimeout(() => setJustSaved(false), 1600);
  };
  const loadSavedDesign = (d) => {
    setWorkspaceOpen(false);
    if (d.tab) setTab(d.tab);
    const r = d.restore || {};
    if (r.seq != null) setSeq(r.seq);
    if (r.seqInfo !== undefined) setSeqInfo(r.seqInfo);
    if (r.mirnaHit !== undefined) setMirnaHit(r.mirnaHit);
    if (r.selectedPair != null) setSelectedPair(r.selectedPair);
    if (r.exonStrategy) setExonStrategy(r.exonStrategy);
    if (r.selectedExons) setSelectedExons(r.selectedExons);
    if (r.method) setMethod(r.method);
    setResult(null); // user re-runs Design with restored inputs
  };

  return (
    <div className="min-h-screen bg-slate-50 text-slate-900">
      <WorkspaceDrawer
        open={workspaceOpen}
        onClose={() => setWorkspaceOpen(false)}
        designs={savedDesigns}
        onLoad={loadSavedDesign}
        onRemove={removeDesign}
        onClear={clearDesigns}
        onExportAll={() => savedDesigns.length && exportPanelCSV(savedDesigns)}
      />
      <div className="mx-auto max-w-5xl px-4 py-8 sm:py-10">

        {/* Header */}
        <header className="mb-7">
          <div className="flex items-start justify-between gap-3">
            <div>
              <div className="flex items-center gap-2 text-teal-700">
                <Dna size={18} />
                <span className="text-[11px] font-semibold uppercase tracking-[0.18em]">FreePrimers</span>
              </div>
              <h1 className="mt-1 text-2xl font-bold tracking-tight text-slate-900 sm:text-3xl">
                Free. Open. Assay-aware primer design.
              </h1>
            </div>
            <div className="mt-1 flex shrink-0 items-center gap-2">
              <button onClick={() => setWorkspaceOpen(true)}
                className="relative inline-flex items-center gap-1.5 rounded-full bg-white px-2.5 py-1 text-[11px] font-medium text-slate-600 ring-1 ring-inset ring-slate-200 hover:bg-slate-50">
                <Bookmark size={12} /> Saved
                {savedDesigns.length > 0 && (
                  <span className="ml-0.5 rounded-full bg-teal-600 px-1.5 text-[10px] font-bold text-white">{savedDesigns.length}</span>
                )}
              </button>
              <button onClick={() => setShowBackendPanel(!showBackendPanel)}
                className={`inline-flex items-center gap-1.5 rounded-full px-2.5 py-1 text-[11px] font-medium ring-1 ring-inset transition-colors ${backendActive ? "bg-teal-50 text-teal-700 ring-teal-600/20" : "bg-white text-slate-500 ring-slate-200 hover:bg-slate-50"}`}>
                {backendActive ? <PlugZap size={12} /> : <Plug size={12} />}
                {backendActive ? "Backend connected" : "Connect backend"}
              </button>
            </div>
          </div>
          <p className="mt-1.5 text-sm text-slate-500">
            Paste a sequence, get ranked primers. Tm from SantaLucia 1998 nearest-neighbor thermodynamics with salt and Mg<sup>2+</sup> corrections.{backendActive ? " Backend connected — values upgrade to real primer3-py thermodynamics automatically." : ""}
          </p>

          {showBackendPanel && (
            <div className="mt-3 rounded-lg border border-slate-200 bg-white p-4">
              <div className="flex items-center justify-between">
                <span className="text-xs font-semibold text-slate-700">FreePrimers backend</span>
                {backendActive && (
                  <button onClick={disconnectBackend} className="text-[11px] font-medium text-slate-400 hover:text-slate-600">Disconnect</button>
                )}
              </div>
              <p className="mt-1 text-xs text-slate-500">
                Optional. Connect a running FreePrimers backend (FastAPI + primer3-py + ViennaRNA + BLAST) to upgrade Tm, hairpin, and dimer values to real nearest-neighbour thermodynamics, and to enable specificity checking. Without one, the tool works exactly as it does now using built-in estimates.
              </p>
              <div className="mt-2 flex gap-2">
                <input value={backendUrlInput} onChange={(e) => setBackendUrlInput(e.target.value)}
                  placeholder="http://localhost:8000" spellCheck={false}
                  className="flex-1 rounded-md border border-slate-200 bg-slate-50 px-2.5 py-1.5 text-sm font-mono outline-none focus:border-teal-500 focus:ring-2 focus:ring-teal-500/20" />
                <button onClick={connectBackend} disabled={backendStatus === "checking"}
                  className="inline-flex items-center gap-1.5 rounded-md bg-slate-900 px-3 py-1.5 text-xs font-semibold text-white hover:bg-slate-800 disabled:bg-slate-300">
                  {backendStatus === "checking" ? <Loader2 size={12} className="animate-spin" /> : <Plug size={12} />}
                  Connect
                </button>
              </div>
              {backendStatus === "connected" && (
                <p className="mt-2 inline-flex items-center gap-1 text-[11px] font-medium text-teal-700"><ShieldCheck size={12} /> Connected. Values will upgrade automatically as designs load.</p>
              )}
              {backendStatus === "failed" && (
                <p className="mt-2 inline-flex items-center gap-1 text-[11px] font-medium text-rose-600"><ShieldAlert size={12} /> Couldn't reach that URL's /health endpoint. Check the backend is running and the URL is correct.</p>
              )}
            </div>
          )}
        </header>

        {/* Tabs */}
        <div className="mb-4 inline-flex rounded-lg border border-slate-200 bg-white p-1">
          {[["mirna", "miRNA", Dna], ["standard", "Standard qPCR / PCR", FlaskConical], ["snp", "SNP genotyping", Beaker], ["placement", "SNP-aware placement", Sparkles], ["hrm", "HRM", FlaskConical], ["methyl", "Methylation", Beaker]].map(([id, lab, Ic]) => (
            <button key={id} onClick={() => reset(id)}
              className={`inline-flex items-center gap-1.5 rounded-md px-3 py-1.5 text-sm font-medium transition-colors ${tab === id ? "bg-slate-900 text-white" : "text-slate-600 hover:bg-slate-100"}`}>
              <Ic size={15} /> {lab}
            </button>
          ))}
        </div>

        {/* Input panel */}
        {tab === "placement" ? (
          <SnpAwarePanel backendActive={backendActive} />
        ) : tab === "hrm" ? (
          <HrmPanel backendActive={backendActive} onSaveDesign={saveExternalDesign} justSaved={justSaved} />
        ) : tab === "methyl" ? (
          <MethylationPanel backendActive={backendActive} onSaveDesign={saveExternalDesign} justSaved={justSaved} />
        ) : (
        <div className="rounded-xl border border-slate-200 bg-white p-5 shadow-sm">
          {tab === "mirna" && (
            <div className="mb-4">
              <div className="mb-2 flex flex-wrap items-center gap-2">
                <span className="text-xs font-medium text-slate-500">Method</span>
                {[["stemloop", "Stem-loop RT-qPCR"], ["polya", "Poly(A) + universal reverse"]].map(([id, lab]) => (
                  <button key={id} onClick={() => setMethod(id)}
                    className={`rounded-full px-3 py-1 text-xs font-medium ring-1 ring-inset transition-colors ${method === id ? "bg-teal-600 text-white ring-teal-600" : "bg-white text-slate-600 ring-slate-200 hover:bg-slate-50"}`}>
                    {lab}
                  </button>
                ))}
              </div>
            </div>
          )}

          {tab === "standard" && (
            <div className="mb-2">
              <CdnaFetchRow setSeq={setSeq} backendActive={backendActive}
                onInfo={(info) => { setSeqInfo(info); setMirnaHit(null); setSelectedExons([]); }} />
              <SequenceInfoCard info={seqInfo}
                exonStrategy={exonStrategy}
                onExonStrategyChange={setExonStrategy}
                selectedExons={selectedExons}
                onSelectExon={((ex) => {
                  if (ex === null) { setSelectedExons([]); return; }
                  setSelectedExons(prev =>
                    prev.includes(ex.exon_number)
                      ? prev.filter(n => n !== ex.exon_number)
                      : [...prev, ex.exon_number].slice(-2)
                  );
                })} />
            </div>
          )}

          {tab === "snp" && (
            <div className="mb-2">
              <SeqLoadRow seq={seq} setSeq={setSeq} backendActive={backendActive}
                lookupMode="rsid"
                onInfo={(info) => { setSeqInfo(info); setMirnaHit(null); }} />
              <SequenceInfoCard info={seqInfo} />
            </div>
          )}

          <div className="relative">
            <textarea
              value={seq} onChange={(e) => setSeq(e.target.value)} rows={tab === "mirna" ? 2 : 4}
              placeholder={placeholder} spellCheck={false}
              className="w-full resize-y rounded-lg border border-slate-200 bg-slate-50 p-3 pr-9 font-mono text-sm text-slate-800 outline-none focus:border-teal-500 focus:ring-2 focus:ring-teal-500/20"
            />
            {seq.length > 0 && (
              <button onClick={() => { setSeq(""); setResult(null); }} title="Clear sequence"
                className="absolute right-2 top-2 rounded-md p-1 text-slate-400 hover:bg-slate-200 hover:text-slate-700">
                <X size={14} />
              </button>
            )}
          </div>

          {tab === "mirna" && (
            <MirnaDbSearch onSelect={(seq, hit) => { setSeq(seq); setMirnaHit(hit); setSeqInfo(null); }} />
          )}
          {tab === "mirna" && <MirnaInfoCard hit={mirnaHit} />}

          {tab === "snp" && (
            <div className="mt-2">
              <div className="mb-2 flex flex-wrap items-center gap-2">
                {/* FASTA upload */}
                <label className="inline-flex cursor-pointer items-center gap-1.5 rounded-md bg-slate-100 px-2.5 py-1 text-[11px] font-medium text-slate-600 hover:bg-slate-200">
                  <Upload size={12} /> Upload FASTA
                  <input type="file" accept=".fa,.fasta,.txt,.seq" className="hidden"
                    onChange={(e) => {
                      const f = e.target.files?.[0];
                      if (!f) return;
                      const reader = new FileReader();
                      reader.onload = () => {
                        const text = String(reader.result || "");
                        const lines = text.split(/\r?\n/);
                        let s = "";
                        for (const ln of lines) {
                          if (ln.startsWith(">")) { if (s) break; continue; }
                          s += ln.trim();
                        }
                        s = s.toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, "");
                        setSeq(s);
                      };
                      reader.readAsText(f);
                      e.target.value = "";
                    }} />
                </label>

                {/* rsID lookup now handled by SeqLoadRow above */}
              </div>
              {/* failure note */}
              {rsidStatus === "failed" && rsidNote && (
                <p className="mb-2 text-[11px] text-amber-600">{rsidNote}</p>
              )}
              {!backendActive && (
                <p className="mb-2 text-[11px] text-slate-400">rsID auto-fetch needs a connected backend (it queries Ensembl). Without one, upload a FASTA or paste the bracketed sequence.</p>
              )}
              <button onClick={() => setSeq(SNP_EXAMPLE)}
                className="rounded-md bg-slate-100 px-2 py-1 text-[11px] font-medium text-slate-600 hover:bg-slate-200">
                Load example SNP region
              </button>

              {/* Click-to-mark SNP picker */}
              {(() => {
                const clean = seq.toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, "");
                const hasBracket = /\[[ACGT]\/[ACGT]\]/i.test(seq);
                if (clean.length < 20) return null;
                if (hasBracket) {
                  return (
                    <p className="mt-2 text-[11px] text-teal-700">
                      SNP marked. <button onClick={() => setSeq(clean)} className="underline hover:text-teal-900">Clear marker</button> to pick a different position.
                    </p>
                  );
                }
                const markSnp = (i) => {
                  const ref = clean[i];
                  const alt = window.prompt(`Mark position ${i + 1} (base ${ref}) as the SNP.\nEnter the second allele (A/C/G/T):`, ref === "C" ? "T" : ref === "A" ? "G" : "A");
                  if (!alt) return;
                  const a = alt.trim().toUpperCase();
                  if (!"ACGT".includes(a) || a.length !== 1) { window.alert("Second allele must be a single base A/C/G/T."); return; }
                  setSeq(clean.slice(0, i) + `[${ref}/${a}]` + clean.slice(i + 1));
                };
                return (
                  <div className="mt-2">
                    <p className="mb-1 text-[11px] text-slate-500">Or click the variant base below to mark it as the SNP:</p>
                    <div className="max-h-28 overflow-y-auto rounded-md border border-slate-200 bg-slate-50 p-2 font-mono text-[11px] leading-5 tracking-wide">
                      {clean.split("").map((b, i) => (
                        <span key={i} onClick={() => markSnp(i)}
                          title={`position ${i + 1}`}
                          className="cursor-pointer rounded-sm px-[0.5px] hover:bg-teal-200 hover:text-teal-900">
                          {b}
                        </span>
                      ))}
                    </div>
                  </div>
                );
              })()}
            </div>
          )}

          {/* Advanced */}
          <button onClick={() => setAdv(!adv)}
            className="mt-3 inline-flex items-center gap-1.5 text-xs font-medium text-slate-500 hover:text-slate-800">
            <Settings2 size={13} /> Advanced settings
            <ChevronDown size={13} className={`transition-transform ${adv ? "rotate-180" : ""}`} />
          </button>
          {adv && (
            <div className="mt-3 grid grid-cols-2 gap-x-5 gap-y-3 rounded-lg bg-slate-50 p-4 sm:grid-cols-3">
              {[
                ["target", "Target Tm (°C)", 1, 50, 72],
                ["na", "[Na⁺] (mM)", 1, 0, 200],
                ["mg", "[Mg²⁺] (mM)", 0.1, 0, 10],
                ["dntp", "[dNTP] (mM)", 0.1, 0, 5],
                ["primer", "[Primer] (nM)", 10, 50, 1000],
              ].map(([k, lab, st, mn, mx]) => (
                <label key={k} className="flex flex-col gap-1">
                  <span className="text-[11px] font-medium text-slate-500">{lab}</span>
                  <input type="number" step={st} min={mn} max={mx} value={p[k]}
                    onChange={(e) => setP({ ...p, [k]: parseFloat(e.target.value) || 0 })}
                    className="rounded-md border border-slate-200 bg-white px-2 py-1 font-mono text-sm outline-none focus:border-teal-500" />
                </label>
              ))}
              {tab === "standard" && (
                <label className="flex flex-col gap-1 col-span-2 sm:col-span-1">
                  <span className="text-[11px] font-medium text-slate-500">Amplicon (bp)</span>
                  <div className="flex items-center gap-1">
                    <input type="number" value={p.amplicon[0]} onChange={(e) => setP({ ...p, amplicon: [parseInt(e.target.value) || 0, p.amplicon[1]] })}
                      className="w-full rounded-md border border-slate-200 bg-white px-2 py-1 font-mono text-sm outline-none focus:border-teal-500" />
                    <span className="text-slate-400">–</span>
                    <input type="number" value={p.amplicon[1]} onChange={(e) => setP({ ...p, amplicon: [p.amplicon[0], parseInt(e.target.value) || 0] })}
                      className="w-full rounded-md border border-slate-200 bg-white px-2 py-1 font-mono text-sm outline-none focus:border-teal-500" />
                  </div>
                </label>
              )}
            </div>
          )}

          <button onClick={run} disabled={!cleanDNA(seq).length}
            className="mt-4 inline-flex w-full items-center justify-center gap-2 rounded-lg bg-slate-900 px-4 py-2.5 text-sm font-semibold text-white transition-colors hover:bg-slate-800 disabled:cursor-not-allowed disabled:bg-slate-300 sm:w-auto">
            <Sparkles size={15} /> Design primers
          </button>
        </div>
        )}

        {/* Results */}
        {result && (
          <div className="mt-6">
            {result.error ? (
              <div className="flex items-start gap-2 rounded-lg border border-amber-200 bg-amber-50 p-4 text-sm text-amber-800">
                <TriangleAlert size={16} className="mt-0.5 shrink-0" /> {result.error}
              </div>
            ) : (
              <>
                <div className="mb-3 flex items-center gap-2">
                  <Beaker size={16} className="text-teal-700" />
                  <h2 className="text-sm font-semibold text-slate-800">{result.method}</h2>
                </div>
                {result.blurb && <p className="mb-3 text-xs leading-relaxed text-slate-500">{result.blurb}</p>}

                {result.oligos && (
                  <>
                    <div className="space-y-3">
                      {result.oligos.map((o, i) => (
                        <OligoCard key={i} o={o} p={result.p} hoverKey={hoverKey} onHover={setHoverKey} backendActive={backendActive} />
                      ))}
                    </div>
                    <div className="mt-3 flex flex-wrap gap-3 text-[11px] text-slate-500">
                      {(result.arms
                        ? [["mir", "primer body"], ["overhang", "deliberate mismatch"], ["tail", "allele-specific 3′ base"], ["scaffold", "control primer"]]
                        : Object.entries(SEG_LABEL)
                      ).map(([k, v]) => (
                        <span key={k} className="inline-flex items-center gap-1">
                          <span className={`inline-block h-3 w-3 rounded ${SEG[k]}`} /> {v}
                        </span>
                      ))}
                    </div>
                  </>
                )}

                {result.bands && (
                  <div className="mt-4 rounded-lg border border-slate-200 bg-white p-4">
                    <div className="mb-2 text-xs font-semibold text-slate-700">Expected gel bands</div>
                    {result.bands.map((b, i) => (
                      <div key={i} className="flex items-center justify-between border-b border-slate-100 py-1.5 text-sm last:border-0">
                        <span className="text-slate-600">{b.name}</span>
                        <span className="font-mono text-slate-900">{b.size} bp</span>
                      </div>
                    ))}
                  </div>
                )}

                {result.pairs && (
                  <div className="space-y-3">
                    {result.pairs.map((pr, i) => (
                      <PairCard key={i} pair={pr} idx={i} active={i === selectedPair} onSelect={() => setSelectedPair(i)} hoverKey={hoverKey} onHover={setHoverKey} backendActive={backendActive} salts={p} exons={seqInfo?.exons} />
                    ))}
                  </div>
                )}

                {(() => {
                  const activeMap = result.pairs ? result.pairs[selectedPair]?.map : result.map;
                  if (!activeMap) return null;
                  const sites = findRestrictionSites(activeMap.seqStr, activeEnzymes);
                  return (
                    <div className="mt-4 rounded-lg border border-slate-200 bg-white p-4">
                      <div className="mb-2 flex items-center justify-between">
                        <span className="text-xs font-semibold text-slate-700">Sequence map</span>
                        {result.pairs && <span className="text-[11px] text-slate-400">Pair {selectedPair + 1}</span>}
                      </div>
                      <AnnotationMap map={activeMap} id="primer-map-svg" hoverKey={hoverKey} onHoverFeature={setHoverKey} sites={sites} />
                      <PrimerAlignment seqStr={activeMap.seqStr} features={activeMap.features} snp={activeMap.snp} exons={seqInfo?.exons} isCdna={seqInfo?.type === "cdna"} />
                      <div className="mt-3 flex flex-wrap items-center gap-1.5 border-t border-slate-100 pt-3">
                        <span className="mr-1 text-[11px] font-medium text-slate-500">Restriction sites</span>
                        {ENZYMES.map((enz) => {
                          const on = activeEnzymes.includes(enz.name);
                          return (
                            <button key={enz.name}
                              onClick={() => setActiveEnzymes(on ? activeEnzymes.filter((n) => n !== enz.name) : [...activeEnzymes, enz.name])}
                              className={`rounded-full px-2 py-0.5 text-[11px] font-medium ring-1 ring-inset transition-colors ${on ? "bg-violet-600 text-white ring-violet-600" : "bg-white text-slate-500 ring-slate-200 hover:bg-slate-50"}`}>
                              {enz.name}
                            </button>
                          );
                        })}
                      </div>
                      {activeMap.seqStr && activeMap.seqStr.length > 4000 && (
                        <p className="mt-2 text-[11px] text-slate-400">Sequence is long — site scan covers the full template.</p>
                      )}
                    </div>
                  );
                })()}

                {result.warn && (
                  <p className="mt-2 text-xs text-amber-700">{result.warn}</p>
                )}

                <SpecificityPanel
                  oligos={getExportOligos(result, selectedPair)}
                  backendActive={backendActive}
                  database={specDatabase}
                  onDatabaseChange={setSpecDatabase}
                  isMirna={tab === "mirna"}
                />

                <ExportBar result={result} mapId="primer-map-svg" selectedPair={selectedPair} onSave={saveCurrentDesign} justSaved={justSaved} />

                <div className="mt-4 flex items-start gap-2 rounded-lg border border-slate-200 bg-slate-100/70 p-3 text-xs text-slate-500">
                  <Info size={14} className="mt-0.5 shrink-0 text-slate-400" />
                  <span>
                    Thermodynamics run locally. <strong className="text-slate-600">Specificity is not checked here</strong> — running these against miRBase and the human transcriptome (BLAST) is the backend step that turns this into a high-accuracy pipeline. Dimer and hairpin values are heuristic estimates, not full minimum-free-energy folding.{result.arms ? " Inner-primer Tm shown is the nominal oligo Tm; allele-specific priming runs slightly cooler because of the deliberate −2 mismatch." : ""}
                  </span>
                </div>
              </>
            )}
          </div>
        )}

        {tab === "snp" && result && !result.error && result.map && (
          <ComparisonPanel
            sense={result.map.seqStr}
            snpIndex={result.map.snp}
            allele1={result.arms?.a1}
            allele2={result.arms?.a2}
            ourDesign={result}
            backendActive={backendActive}
          />
        )}

        {tab === "mirna" && (
          <MirnaComparePanel backendActive={backendActive} />
        )}

        {!result && tab !== "placement" && tab !== "hrm" && tab !== "methyl" && (
          <div className="mt-6 rounded-xl border border-dashed border-slate-200 p-8 text-center">
            <Dna size={22} className="mx-auto text-slate-300" />
            <p className="mt-2 text-sm text-slate-400">
              {tab === "mirna" ? "Paste a mature miRNA or pick an example above to design primers." : "Paste a target sequence to find candidate primer pairs."}
            </p>
          </div>
        )}

        <footer className="mt-8 text-center text-[11px] leading-relaxed text-slate-400">
          FreePrimers — open-source, assay-aware primer design. Always validate empirically. Stem-loop scaffold follows Chen / Varkonyi-Gasic; miRNA presets are examples — confirm against miRBase. MIT License.
        </footer>
      </div>
    </div>
  );
}
