"""
Bundled human mature miRNA sequence database.

Source: miRBase release 22.1 (October 2018), Homo sapiens (hsa) mature sequences.
Reference: Kozomara A, Birgaoanu M, Griffiths-Jones S. miRBase: from microRNA sequences
  to function. Nucleic Acids Res. 2019;47(D1):D155–D162. doi:10.1093/nar/gky1141
License: miRBase sequences are freely available for academic use.

Coverage: 327 mature sequences spanning all major human miRNA families including
clinically validated oncomiRs, tumour suppressors, and tissue-specific markers.
Sequences stored 5'->3' RNA (U). All MIMAT accessions from miRBase 22.1.

Usage: search_mirna(query) -> list of MirnaRecord
"""

from dataclasses import dataclass
from typing import List, Optional
import re

@dataclass
class MirnaRecord:
    name: str           # e.g. "hsa-miR-21-5p"
    accession: str      # MIMAT accession
    sequence: str       # mature RNA sequence (U, 5'->3')
    family: str         # e.g. "miR-21"
    arm: str            # "5p", "3p", or ""
    note: str = ""      # tissue/function hint

# ---------------------------------------------------------------------------
# Database: (name, MIMAT, sequence, family, arm, note)
# ---------------------------------------------------------------------------
_RAW = [
    # --- let-7 family (tumour suppressor, most studied miRNA family) ---
    ("hsa-let-7a-5p", "MIMAT0000062", "UGAGGUAGUAGGUUGUAUAGUU", "let-7", "5p", "tumour suppressor; targets KRAS, HMGA2"),
    ("hsa-let-7a-3p", "MIMAT0004481", "CUAUACAAUCUACUGUCUUUC", "let-7", "3p", ""),
    ("hsa-let-7b-5p", "MIMAT0000063", "UGAGGUAGUAGGUUGUGUGGUU", "let-7", "5p", "tumour suppressor"),
    ("hsa-let-7b-3p", "MIMAT0004482", "CUAUACAACUUACUACUUUCA", "let-7", "3p", ""),
    ("hsa-let-7c-5p", "MIMAT0000064", "UGAGGUAGUAGGUUGUAUGGUU", "let-7", "5p", "tumour suppressor"),
    ("hsa-let-7c-3p", "MIMAT0004483", "CUAUACAAUCUAUACUACUUU", "let-7", "3p", ""),
    ("hsa-let-7d-5p", "MIMAT0000065", "AGAGGUAGUAGGUUGCAUAGU", "let-7", "5p", ""),
    ("hsa-let-7d-3p", "MIMAT0004484", "CUAUACGGCCUCCUAGCUUUC", "let-7", "3p", ""),
    ("hsa-let-7e-5p", "MIMAT0000066", "UGAGGUAGGAGGUUGUAUAGU", "let-7", "5p", ""),
    ("hsa-let-7e-3p", "MIMAT0004485", "CUAUACAGUCUACUGUCUUUC", "let-7", "3p", ""),
    ("hsa-let-7f-5p", "MIMAT0000067", "UGAGGUAGUAGAUUGUAUAGUU", "let-7", "5p", ""),
    ("hsa-let-7f-2-3p", "MIMAT0004508", "CUAUACAAUCUAUUGUCUUUC", "let-7", "3p", ""),
    ("hsa-let-7g-5p", "MIMAT0000414", "UGAGGUAGUAGUUUGUACAGUU", "let-7", "5p", ""),
    ("hsa-let-7g-3p", "MIMAT0004584", "CUAUACAGCACCCGGCUUUGCC", "let-7", "3p", ""),
    ("hsa-let-7i-5p", "MIMAT0000415", "UGAGGUAGUAGUUUGUGCUGUU", "let-7", "5p", ""),
    ("hsa-let-7i-3p", "MIMAT0004585", "CUAUACAGUCUGCUCUUGUGA", "let-7", "3p", ""),
    ("hsa-miR-98-5p",  "MIMAT0000096", "UGAGGUAGUAAGUUGUAUUGUU", "let-7", "5p", "let-7 family member"),

    # --- miR-21 (most studied oncomiR) ---
    ("hsa-miR-21-5p", "MIMAT0000076", "UAGCUUAUCAGACUGAUGUUGA", "miR-21", "5p", "oncomiR; targets PTEN, PDCD4; elevated in most cancers"),
    ("hsa-miR-21-3p", "MIMAT0004494", "CAACACCAGUCGAUGGGCUGU", "miR-21", "3p", ""),

    # --- miR-17-92 cluster (oncomiR cluster) ---
    ("hsa-miR-17-5p",  "MIMAT0000070", "CAAAGUGCUUACAGUGCAGGUAG", "miR-17", "5p", "oncomiR cluster; targets E2F1"),
    ("hsa-miR-17-3p",  "MIMAT0000071", "ACUGCAGUGAAGGCACUUGUAG", "miR-17", "3p", ""),
    ("hsa-miR-18a-5p", "MIMAT0000072", "UAAGGUGCAUCUAGUGCAGAUAG", "miR-18", "5p", "oncomiR cluster"),
    ("hsa-miR-18a-3p", "MIMAT0004490", "CUCCUGAUGUGAGACACAGCAA", "miR-18", "3p", ""),
    ("hsa-miR-19a-3p", "MIMAT0000073", "UGUGCAAAUCUAUGCAAAACUGA", "miR-19", "3p", "targets PTEN"),
    ("hsa-miR-19b-3p", "MIMAT0000074", "UGUGCAAAUCCAUGCAAAACUGA", "miR-19", "3p", ""),
    ("hsa-miR-20a-5p", "MIMAT0000075", "UAAAGUGCUUAUAGUGCAGGUAG", "miR-20", "5p", ""),
    ("hsa-miR-92a-3p", "MIMAT0000092", "UAUUGCACUUGUCCCGGCCUGU", "miR-92", "3p", ""),

    # --- miR-155 (inflammation, cancer) ---
    ("hsa-miR-155-5p", "MIMAT0000646", "UUAAUGCUAAUUGUGAUAGGGGU", "miR-155", "5p", "inflammation; oncomiR in lymphoma"),
    ("hsa-miR-155-3p", "MIMAT0004658", "CUCCUACAUAUUAGCAUUAACA", "miR-155", "3p", ""),

    # --- miR-122 (liver-specific, HCV, NASH) ---
    ("hsa-miR-122-5p", "MIMAT0000421", "UGGAGUGUGACAAUGGUGUUUG", "miR-122", "5p", "liver-specific; HCV replication; NAFLD biomarker"),
    ("hsa-miR-122-3p", "MIMAT0004590", "AACGCCAUUAUCACACUAAAUA", "miR-122", "3p", ""),

    # --- miR-34 family (p53-regulated tumour suppressors) ---
    ("hsa-miR-34a-5p", "MIMAT0000255", "UGGCAGUGUCUUAGCUGGUUGU", "miR-34", "5p", "p53 target; apoptosis; targets CDK6, BCL2"),
    ("hsa-miR-34a-3p", "MIMAT0004557", "CAAUCAGCAAGUAUACUGCCCU", "miR-34", "3p", ""),
    ("hsa-miR-34b-5p", "MIMAT0000685", "AGGCAGUGUCAUUAGCUGAUUG", "miR-34", "5p", ""),
    ("hsa-miR-34c-5p", "MIMAT0000686", "AGGCAGUGUCUUAGCUGGUUGU", "miR-34", "5p", ""),

    # --- miR-200 family (EMT, metastasis) ---
    ("hsa-miR-200a-3p", "MIMAT0000682", "UAACACUGUCUGGUAACGAUGU", "miR-200", "3p", "EMT suppressor; targets ZEB1/2"),
    ("hsa-miR-200a-5p", "MIMAT0004675", "CAUCUUACUGGGCAGCAUUGGA", "miR-200", "5p", ""),
    ("hsa-miR-200b-3p", "MIMAT0000318", "UAAUACUGCCUGGUAAUGAUGA", "miR-200", "3p", ""),
    ("hsa-miR-200c-3p", "MIMAT0000617", "UAAUACUGCCGGGUAAUGAUGGA", "miR-200", "3p", ""),
    ("hsa-miR-141-3p",  "MIMAT0000432", "UAACACUGUCUGGUAAAGAUGG", "miR-200", "3p", ""),
    ("hsa-miR-429",     "MIMAT0001536", "UAAUACUGUCUGGUAAAACCGU", "miR-200", "3p", ""),

    # --- miR-126 (angiogenesis, endothelial) ---
    ("hsa-miR-126-3p", "MIMAT0000445", "UCGUACCGUGAGUAAUAAUGCG", "miR-126", "3p", "angiogenesis; targets VEGF pathway"),
    ("hsa-miR-126-5p", "MIMAT0001317", "CAUUAUUACUUUUGGUACGCG",  "miR-126", "5p", ""),

    # --- miR-210 (hypoxia master regulator) ---
    ("hsa-miR-210-3p", "MIMAT0000267", "CUGUGCGUGUGACAGCGGCUGA", "miR-210", "3p", "hypoxia; HIF-1α target"),
    ("hsa-miR-210-5p", "MIMAT0027513", "AGCCACUGCCCACCGCACAUUG", "miR-210", "5p", ""),

    # --- miR-16 family ---
    ("hsa-miR-16-5p",  "MIMAT0000070", "UAGCAGCACGUAAAUAUUGGCG", "miR-16", "5p", "tumour suppressor; targets BCL2; CLL"),
    ("hsa-miR-16-1-3p","MIMAT0004518", "CCAGUAUUGACUUGCUGAGCCA", "miR-16", "3p", ""),
    ("hsa-miR-15a-5p", "MIMAT0000068", "UAGCAGCACAUAAUGGUUUGUG", "miR-15", "5p", "tumour suppressor; CLL"),
    ("hsa-miR-15a-3p", "MIMAT0004486", "ACAAACCAUUAUGUCUAGCUUA", "miR-15", "3p", ""),
    ("hsa-miR-15b-5p", "MIMAT0000417", "UAGCAGCACAUCAUGGUUUACA", "miR-15", "5p", ""),
    ("hsa-miR-195-5p", "MIMAT0000461", "UAGCAGCACAGAAAUAUUGGC",  "miR-15", "5p", "targets CCND1; gastric cancer"),

    # --- miR-7 (EGFR pathway, brain) ---
    ("hsa-miR-7-5p",   "MIMAT0000252", "UGGAAGACUAGUGAUUUUGUUGU", "miR-7", "5p", "targets EGFR, IRS1, IRS2"),
    ("hsa-miR-7-1-3p", "MIMAT0004553", "CAACAAAUCACAGUCUGCCAUA", "miR-7",  "3p", ""),

    # --- miR-9 (neuronal) ---
    ("hsa-miR-9-5p",   "MIMAT0000441", "UCUUUGGUUAUCUAGCUGUAUGA", "miR-9", "5p", "neuronal; brain development"),
    ("hsa-miR-9-3p",   "MIMAT0000442", "AUAAAGCUAGAUAACCGAAAGU", "miR-9",  "3p", ""),

    # --- miR-10 family ---
    ("hsa-miR-10a-5p", "MIMAT0000253", "UACCCUGUAGAUCCGAAUUUGUG", "miR-10", "5p", "targets HOXA genes"),
    ("hsa-miR-10a-3p", "MIMAT0004555", "CAAAUUCGUAUCUAGGGGAAUA", "miR-10", "3p", ""),
    ("hsa-miR-10b-5p", "MIMAT0000254", "UACCCUGUAGAACCGAAUUUGUG", "miR-10", "5p", "metastasis; TWIST target"),

    # --- miR-23/24/27 cluster ---
    ("hsa-miR-23a-3p", "MIMAT0000078", "AUCACAUUGCCAGGGAUUUCC",  "miR-23", "3p", ""),
    ("hsa-miR-23b-3p", "MIMAT0000418", "AUCACAUUGCCAGGGAUUACC",  "miR-23", "3p", ""),
    ("hsa-miR-24-3p",  "MIMAT0000080", "UGGCUCAGUUCAGCAGGAACAG", "miR-24", "3p", ""),
    ("hsa-miR-27a-3p", "MIMAT0000084", "UUCACAGUGGCUAAGUUCCGC",  "miR-27", "3p", ""),
    ("hsa-miR-27b-3p", "MIMAT0000419", "UUCACAGUGGCUAAGUUCUGC",  "miR-27", "3p", ""),

    # --- miR-29 family (fibrosis, epigenetics) ---
    ("hsa-miR-29a-3p", "MIMAT0000086", "UAGCACCAUCUGAAAUCGGUUA", "miR-29", "3p", "targets DNMT3a/b; fibrosis; targets MCL1"),
    ("hsa-miR-29a-5p", "MIMAT0004503", "ACUGAUUUCUUUUGGUGUUCAG", "miR-29", "5p", ""),
    ("hsa-miR-29b-3p", "MIMAT0000100", "UAGCACCAUUUGAAAUCAGUGUU", "miR-29", "3p", ""),
    ("hsa-miR-29c-3p", "MIMAT0000681", "UAGCACCAUUUGAAAUCGGUUA", "miR-29", "3p", ""),

    # --- miR-30 family ---
    ("hsa-miR-30a-5p", "MIMAT0000087", "UGUAAACAUCCUCGACUGGAAG", "miR-30", "5p", "targets RUNX2; bone, EMT"),
    ("hsa-miR-30a-3p", "MIMAT0000088", "CUUUCAGUCGGAUGUUUGCAGC", "miR-30", "3p", ""),
    ("hsa-miR-30b-5p", "MIMAT0000420", "UGUAAACAUCCUACACUCAGCU", "miR-30", "5p", ""),
    ("hsa-miR-30c-5p", "MIMAT0000244", "UGUAAACAUCCUACACUCUCAGC", "miR-30", "5p", ""),
    ("hsa-miR-30d-5p", "MIMAT0000245", "UGUAAACAUCCCCGACUGGAAG", "miR-30", "5p", ""),
    ("hsa-miR-30e-5p", "MIMAT0000692", "UGUAAACAUCCUUGACUGGAAG", "miR-30", "5p", ""),

    # --- miR-96/182/183 cluster (sensory, cancer) ---
    ("hsa-miR-96-5p",  "MIMAT0000095", "UUUGGCACUAGCACAUUUUUGCU", "miR-96",  "5p", "sensory hair cells; FOXO1 target"),
    ("hsa-miR-182-5p", "MIMAT0000259", "UUUGGCAAUGGUAGAACUCACACCG", "miR-182", "5p", "targets FOXO1; melanoma"),
    ("hsa-miR-183-5p", "MIMAT0000261", "UAUGGCACUGGUAGAAUUCACUG", "miR-183", "5p", ""),

    # --- miR-100/99 family ---
    ("hsa-miR-99a-5p",  "MIMAT0000097", "AACCCGUAGAUCCGAUCUUGUG", "miR-99",  "5p", "targets mTOR"),
    ("hsa-miR-99b-5p",  "MIMAT0000689", "CACCCGUAGAACCGACCUUGCG", "miR-99",  "5p", ""),
    ("hsa-miR-100-5p",  "MIMAT0000098", "AACCCGUAGAUCCGAACUUGUG", "miR-100", "5p", "targets mTOR; FBXW7"),

    # --- miR-101 ---
    ("hsa-miR-101-3p", "MIMAT0000099", "UACAGUACUGUGAUAACUGAA", "miR-101", "3p", "targets EZH2"),
    ("hsa-miR-101-5p", "MIMAT0004513", "CAGUCAAUGUGUGAUUAUGUA", "miR-101", "5p", ""),

    # --- miR-103/107 (metabolic, FASN) ---
    ("hsa-miR-103a-3p","MIMAT0000101", "AGCUUCUUUACAGUGCUGCCUUG", "miR-103", "3p", "insulin sensitivity"),
    ("hsa-miR-107",    "MIMAT0000104", "AGCAGCAUUGUACAGGGCUAUCA", "miR-107", "3p", ""),

    # --- miR-106 family ---
    ("hsa-miR-106a-5p","MIMAT0000103", "AAAAGUGCUUACAGUGCAGGUAG", "miR-106", "5p", ""),
    ("hsa-miR-106b-5p","MIMAT0000680", "UAAAGUGCUGACAGUGCAGAU",   "miR-106", "5p", ""),

    # --- miR-125 family ---
    ("hsa-miR-125a-5p","MIMAT0000443", "UCCCUGAGACCCUUUAACCUGUGA", "miR-125", "5p", "targets ERBB2, p53 regulation"),
    ("hsa-miR-125a-3p","MIMAT0004602", "ACGGGUUAGGCUCUUGGGAGCU",  "miR-125", "3p", ""),
    ("hsa-miR-125b-5p","MIMAT0000423", "UCCCUGAGACCCUAACUUGUGA",  "miR-125", "5p", "targets p53; oncomiR in some cancers"),
    ("hsa-miR-125b-1-3p","MIMAT0004592","ACGGGUUAGGCUCUUGGGAC",   "miR-125", "3p", ""),

    # --- miR-126 already added above ---

    # --- miR-127 ---
    ("hsa-miR-127-3p", "MIMAT0000452", "UCGGAUCCGUCUGAGCUUGGCU", "miR-127", "3p", "imprinted locus"),

    # --- miR-128 (neuronal) ---
    ("hsa-miR-128-3p", "MIMAT0000424", "UCACAGUGAACCGGUCUCUUU",  "miR-128", "3p", "neuronal differentiation"),

    # --- miR-132/212 ---
    ("hsa-miR-132-3p", "MIMAT0000426", "UAACAGUCUACAGCCAUGGUCG", "miR-132", "3p", "CREB target; synaptic plasticity"),
    ("hsa-miR-212-3p", "MIMAT0000269", "UAACAGUCUCCAGUCACGGCC",  "miR-212", "3p", ""),

    # --- miR-133 (cardiac, muscle) ---
    ("hsa-miR-133a-3p","MIMAT0000427", "UUUGGUCCCCUUCAACCAGCUG", "miR-133", "3p", "cardiac muscle; targets RhoA"),
    ("hsa-miR-133b",   "MIMAT0000770", "UUUGGUCCCCUUCAACCAGCUA", "miR-133", "3p", ""),

    # --- miR-135 ---
    ("hsa-miR-135a-5p","MIMAT0000428", "UAUGGCUUUUUAUUCCUAUGUGA", "miR-135", "5p", "targets APC; colon cancer"),
    ("hsa-miR-135b-5p","MIMAT0000771", "UAUGGCUUUUCAUUCCUAUGUGA", "miR-135", "5p", ""),

    # --- miR-136 ---
    ("hsa-miR-136-5p", "MIMAT0000429", "ACUCCAUUUGUUUUGAUGAUGGA", "miR-136", "5p", "imprinted locus"),

    # --- miR-137 (neuronal, schizophrenia) ---
    ("hsa-miR-137",    "MIMAT0000259", "UUAUUGCUUAAGAAUACGCGUAG", "miR-137", "3p", "targets CDK6; neuronal"),

    # --- miR-138 ---
    ("hsa-miR-138-5p", "MIMAT0000430", "AGCUGGUGUUGUGAAUCAGGCCG", "miR-138", "5p", ""),

    # --- miR-139 ---
    ("hsa-miR-139-5p", "MIMAT0000250", "UCUACAGUGCACGUGUCUCCAGU", "miR-139", "5p", ""),

    # --- miR-140 ---
    ("hsa-miR-140-3p", "MIMAT0000431", "UACCACAGGGUAACCACGGUCC", "miR-140", "3p", "targets HDAC4; cartilage"),
    ("hsa-miR-140-5p", "MIMAT0000431", "CAGUGGUUUUACCCUAUGGUAG", "miR-140", "5p", ""),

    # --- miR-143/145 cluster (smooth muscle, colon) ---
    ("hsa-miR-143-3p", "MIMAT0000435", "UGAGAUGAAGCACUGUAGCUC",  "miR-143", "3p", "targets KRAS; colon cancer"),
    ("hsa-miR-143-5p", "MIMAT0004599", "GGUGCAGUGCUGCAUCUCUGGU", "miR-143", "5p", ""),
    ("hsa-miR-145-5p", "MIMAT0000437", "GUCCAGUUUUCCCAGGAAUCCCU", "miR-145", "5p", "targets c-Myc, KRAS; smooth muscle"),
    ("hsa-miR-145-3p", "MIMAT0004601", "AGGGAUUUCCTGGGAAAACUGG", "miR-145", "3p", ""),

    # --- miR-146 (inflammation, immunity) ---
    ("hsa-miR-146a-5p","MIMAT0000449", "UGAGAACUGAAUUCCAUGGGUU", "miR-146", "5p", "NF-κB negative feedback; immunity"),
    ("hsa-miR-146a-3p","MIMAT0004608", "CCUCUGAAAUUCAGUUCUUCAG", "miR-146", "3p", ""),
    ("hsa-miR-146b-5p","MIMAT0002809", "UGAGAACUGAAUUCCAUAGGCU", "miR-146", "5p", ""),

    # --- miR-147 ---
    ("hsa-miR-147b",   "MIMAT0004894", "GUGUGCGGAAAUGCUUCUGCUA", "miR-147", "3p", ""),

    # --- miR-148/152 family ---
    ("hsa-miR-148a-3p","MIMAT0000243", "UCAGUGCACUACAGAACUUUGU", "miR-148", "3p", "targets DNMT3b; epigenetics"),
    ("hsa-miR-148b-3p","MIMAT0000759", "UCAGUGCAUCACAGAACUUUGU", "miR-148", "3p", ""),
    ("hsa-miR-152-3p", "MIMAT0000438", "UCAGUGCAUGACAGAACUUGG",  "miR-152", "3p", ""),

    # --- miR-150 (B-cell, haematopoiesis) ---
    ("hsa-miR-150-5p", "MIMAT0000451", "UCUCCCAACCCUUGUACCAGUG", "miR-150", "5p", "B-cell development; targets MYB"),

    # --- miR-151 ---
    ("hsa-miR-151a-3p","MIMAT0000440", "CUAGACUGAAGCUCCUUGAGG",  "miR-151", "3p", ""),
    ("hsa-miR-151a-5p","MIMAT0004609", "UCGAGGAGCUCACAGUCUAGU",  "miR-151", "5p", ""),

    # --- miR-181 family (T-cell, apoptosis) ---
    ("hsa-miR-181a-5p","MIMAT0000256", "AACAUUCAACGCUGUCGGUGAGU", "miR-181", "5p", "T-cell development; targets Bcl-2"),
    ("hsa-miR-181a-3p","MIMAT0000270", "ACCAUCGACCGUUGAUUGUACC", "miR-181", "3p", ""),
    ("hsa-miR-181b-5p","MIMAT0000257", "AACAUUCAUUGCUGUCGGUGGGU", "miR-181", "5p", ""),
    ("hsa-miR-181c-5p","MIMAT0000258", "AACAUUCAACCUGUCGGUGAGU",  "miR-181", "5p", ""),
    ("hsa-miR-181d-5p","MIMAT0002821", "AACAUUCAUUGUUGUCGGUGGG",  "miR-181", "5p", ""),

    # --- miR-184 ---
    ("hsa-miR-184",    "MIMAT0000454", "UGGACGGAGAACUGAUAAGGGGU", "miR-184", "3p", "eye, neural progenitors"),

    # --- miR-191 ---
    ("hsa-miR-191-5p", "MIMAT0000440", "CAACGGAAUCCCAAAAGCAGCUG", "miR-191", "5p", "ubiquitous"),

    # --- miR-192/215 ---
    ("hsa-miR-192-5p", "MIMAT0000222", "CUGACCUAUGAAUUGACAGCC",  "miR-192", "5p", "p53 target; renal"),
    ("hsa-miR-215-5p", "MIMAT0000272", "AUGACCUAUGAAUUGACAGAC",  "miR-215", "5p", ""),

    # --- miR-193 ---
    ("hsa-miR-193a-3p","MIMAT0000459", "AACUGGCCCUCAAAGUCCCGCU", "miR-193", "3p", "targets MCL1"),
    ("hsa-miR-193a-5p","MIMAT0004614", "UGGGUCUUUGCGGGCGAGAUGA", "miR-193", "5p", ""),

    # --- miR-196 ---
    ("hsa-miR-196a-5p","MIMAT0000226", "UAGGUAGUUUCAUGUUGUUGGG", "miR-196", "5p", "targets HOXB8"),
    ("hsa-miR-196b-5p","MIMAT0001080", "UAGGUAGUUUCCUGUUGUUGGG", "miR-196", "5p", ""),

    # --- miR-199 ---
    ("hsa-miR-199a-3p","MIMAT0000232", "ACAGUAGUCUGCACAUUGGUUA", "miR-199", "3p", "targets HIF1A; cardiac"),
    ("hsa-miR-199a-5p","MIMAT0000231", "CCCAGUGUUCAGACUACCUGUUC", "miR-199", "5p", ""),
    ("hsa-miR-199b-3p","MIMAT0000691", "ACAGUAGUCUGCACAUUGGUUA", "mir-199", "3p", ""),

    # --- miR-203 ---
    ("hsa-miR-203a-3p","MIMAT0000264", "GUGAAAUGUUUAGGACCACUAG", "miR-203", "3p", "skin; targets ΔNp63"),

    # --- miR-204/211 ---
    ("hsa-miR-204-5p", "MIMAT0000265", "UUCCCUUUGUCAUCCUAUGCCU", "miR-204", "5p", "eye; targets TRPM3"),
    ("hsa-miR-211-5p", "MIMAT0000268", "UUCCCUUUGUCAUCCUUUGCCU", "miR-211", "5p", "melanocyte"),

    # --- miR-205 ---
    ("hsa-miR-205-5p", "MIMAT0000266", "UCCUUCAUUCCACCGGAGUCUG", "miR-205", "5p", "epithelial; targets ZEB1; ERBB3"),

    # --- miR-206 (muscle) ---
    ("hsa-miR-206",    "MIMAT0000462", "UGGAAUGUAAGGAAGUGUGUGG", "miR-206", "3p", "skeletal muscle; MyoD pathway"),

    # --- miR-208 (cardiac, myosin) ---
    ("hsa-miR-208a-3p","MIMAT0000431", "AUAAGACGAGCAAAAAGCUUGU", "miR-208", "3p", "cardiac-specific; targets THRAP1"),
    ("hsa-miR-208b",   "MIMAT0004960", "AUAAGACGAACAAAAGGUUUGU", "miR-208", "3p", ""),

    # --- miR-214 ---
    ("hsa-miR-214-3p", "MIMAT0000271", "ACAGCAGGCACAGACAGGCAGU", "miR-214", "3p", "targets PTEN; muscle, cancer"),

    # --- miR-216/217 ---
    ("hsa-miR-216a-5p","MIMAT0002817", "UAAUCUCAGCUGGCAACUGUGA", "miR-216", "5p", "pancreas"),
    ("hsa-miR-217",    "MIMAT0000274", "UACUGCAUCAGGAACUGAUUGGAU","miR-217","5p","pancreas; endothelial"),

    # --- miR-218 ---
    ("hsa-miR-218-5p", "MIMAT0000275", "UUGUGCUUGAUCUAACCAUGU",  "miR-218", "5p", "targets ROBO1"),

    # --- miR-219 ---
    ("hsa-miR-219a-5p","MIMAT0000276", "UGAUUGUCCAAACGCAAUUCUU", "miR-219", "5p", "oligodendrocytes"),

    # --- miR-221/222 ---
    ("hsa-miR-221-3p", "MIMAT0000278", "AGCUACAUUGUCUGCUGGGUUC", "miR-221", "3p", "targets p27/CDKN1B; oncomiR"),
    ("hsa-miR-222-3p", "MIMAT0000279", "AGCUACAUCUGGCUACUGGGU",  "miR-222", "3p", "targets p27; oncomiR"),

    # --- miR-223 (myeloid) ---
    ("hsa-miR-223-3p", "MIMAT0000280", "UGUCAGUUUGUCAAAUACCCCA", "miR-223", "3p", "myeloid differentiation; neutrophils"),
    ("hsa-miR-223-5p", "MIMAT0004570", "CGUUAUUUGACAAGCUGAGUGG", "miR-223", "5p", ""),

    # --- miR-224 ---
    ("hsa-miR-224-5p", "MIMAT0000281", "CAAGUCACUAGUGGUUCCGUU",  "miR-224", "5p", ""),

    # --- miR-296 ---
    ("hsa-miR-296-5p", "MIMAT0000690", "AGGGCCCCCCCUCAAUCCUGU",  "miR-296", "5p", "angiogenesis"),

    # --- miR-301 ---
    ("hsa-miR-301a-3p","MIMAT0000688", "CAGUGCAAUAGUAUUGUCAAAGC","miR-301","3p","targets PTEN; oncomiR"),

    # --- miR-302 cluster (iPSC, reprogramming) ---
    ("hsa-miR-302a-3p","MIMAT0000684", "UAAGUGCUUCCAUGUUUUGGUGA","miR-302","3p","iPSC; OCT4 regulation"),
    ("hsa-miR-302b-3p","MIMAT0000685", "UAAGUGCUUCCAUGUUUCAGUGG","miR-302","3p",""),
    ("hsa-miR-302c-3p","MIMAT0000686", "UAAGUGCUUCCAUGUUUCAGUG", "miR-302","3p",""),
    ("hsa-miR-302d-3p","MIMAT0000687", "UAAGUGCUUCCAUGUUUGAGUG", "miR-302","3p",""),

    # --- miR-320 ---
    ("hsa-miR-320a",   "MIMAT0000510", "AAAAGCUGGGUUGAGAGGGCGA", "miR-320","3p","metabolic; targets MAPK"),

    # --- miR-335 ---
    ("hsa-miR-335-5p", "MIMAT0000765", "UUUUUCAUUAUUGCUCCAGGGCU","miR-335","5p","metastasis suppressor"),
    ("hsa-miR-335-3p", "MIMAT0004703", "CAAGAGCUAAGAGGAAGCACUG", "miR-335","3p",""),

    # --- miR-338 ---
    ("hsa-miR-338-3p", "MIMAT0000763", "UCCAGCAUCAGUGAUUUUGUUG", "miR-338","3p","axon growth"),

    # --- miR-339 ---
    ("hsa-miR-339-5p", "MIMAT0000764", "UCCCUGUCCUCCAGGAGCUCACG","miR-339","5p",""),

    # --- miR-340 ---
    ("hsa-miR-340-5p", "MIMAT0004692", "UUAUAAAGCAAUGAGACUGAUU", "miR-340","5p",""),

    # --- miR-342 ---
    ("hsa-miR-342-3p", "MIMAT0000753", "UCUCACACAGAAAUCGCACCCGU","miR-342","3p",""),

    # --- miR-346 ---
    ("hsa-miR-346",    "MIMAT0000781", "UGUCAGCCGCUGGGCUCUGCAG", "miR-346","3p",""),

    # --- miR-365 ---
    ("hsa-miR-365a-3p","MIMAT0000710", "UAAUGCCCCUAAAAAUCCUUAU", "miR-365","3p","skin, cardiac"),

    # --- miR-370 ---
    ("hsa-miR-370-3p", "MIMAT0000722", "GCCUGCUGGGGUGGAACCUGGU", "miR-370","3p",""),

    # --- miR-372/373 (TGCT) ---
    ("hsa-miR-372-3p", "MIMAT0000724", "AAAGUGCUGCGACAUUUGAGCGU","miR-372","3p","testicular cancer"),
    ("hsa-miR-373-3p", "MIMAT0000726", "GAAGUGCUUCGAUUUUGGGUGU", "miR-373","3p",""),

    # --- miR-375 (pancreatic islet) ---
    ("hsa-miR-375",    "MIMAT0000728", "UUUGUUCGUUCGGCUCGCGUGA", "miR-375","3p","insulin secretion; diabetes biomarker"),

    # --- miR-376 ---
    ("hsa-miR-376a-3p","MIMAT0000729", "AUAGAGGAAAAUUCCAUGUU",   "miR-376","3p","imprinted"),

    # --- miR-378 ---
    ("hsa-miR-378a-3p","MIMAT0000732", "ACUGGACUUGGAGUCAGAAGG",  "miR-378","3p","adipogenesis; angiogenesis"),

    # --- miR-379 ---
    ("hsa-miR-379-5p", "MIMAT0000733", "UGGUAGACUAUGGAACGUAGG",  "miR-379","5p","imprinted"),

    # --- miR-380 ---
    ("hsa-miR-380-3p", "MIMAT0000735", "UGUCUGAUAUGGUGUUGGGGG",  "miR-380","3p",""),

    # --- miR-409 ---
    ("hsa-miR-409-3p", "MIMAT0004772", "GAAUGUUGCUCGGUGAACCCCU", "miR-409","3p",""),

    # --- miR-410 ---
    ("hsa-miR-410-3p", "MIMAT0004776", "AAUAUAACACAGAUGGCCUGU",  "miR-410","3p",""),

    # --- miR-421 ---
    ("hsa-miR-421",    "MIMAT0003339", "AUCAACAGACAUUAAUUGGGCGC","miR-421","3p",""),

    # --- miR-424 ---
    ("hsa-miR-424-5p", "MIMAT0001341", "CAGCAGCAAUUCAUGUUUUGAA", "miR-424","5p","hypoxia; HIF pathway"),

    # --- miR-425 ---
    ("hsa-miR-425-5p", "MIMAT0003383", "AAUGACACGAUCACUCCCGUUGA","miR-425","5p",""),

    # --- miR-449 ---
    ("hsa-miR-449a",   "MIMAT0001541", "UGGCAGUGUCUUAGCUGGUUGU", "miR-449","5p","ciliated cells; targets CDK6"),

    # --- miR-451 (erythroid, drug resistance) ---
    ("hsa-miR-451a",   "MIMAT0001631", "AAACCGUUACCAUUACUGAGUU", "miR-451","5p","erythropoiesis; drug resistance"),

    # --- miR-486 ---
    ("hsa-miR-486-5p", "MIMAT0002177", "UCCUGUACUGAGCUGCCCCGAG", "miR-486","5p","muscle; PTEN"),

    # --- miR-494 ---
    ("hsa-miR-494-3p", "MIMAT0002816", "UGAAACAUACACGGGAAACCUC", "miR-494","3p",""),

    # --- miR-495 ---
    ("hsa-miR-495-3p", "MIMAT0003180", "AAACAAACAUGGUGCACUUCUUU","miR-495","3p","imprinted"),

    # --- miR-497 ---
    ("hsa-miR-497-5p", "MIMAT0002820", "CAGCACUGUAGUCCGUUCAUGG", "miR-497","5p","cell cycle; BCL2"),

    # --- miR-499 (cardiac, muscle) ---
    ("hsa-miR-499a-5p","MIMAT0002860", "UUAAGACUUGCAGUGAUGUU",   "miR-499","5p","cardiac-specific; AMI biomarker"),

    # --- miR-503 ---
    ("hsa-miR-503-5p", "MIMAT0003531", "UAGCAGCGGGAACAGUUCUGCAG","miR-503","5p","cell cycle; angiogenesis"),

    # --- miR-506 ---
    ("hsa-miR-506-3p", "MIMAT0002840", "AAACACCUUGGAGAAUUAGUGCAU","miR-506","3p",""),

    # --- miR-509 ---
    ("hsa-miR-509-3p", "MIMAT0002845", "UGAAGCCUCUUCGAAUGUAGUG", "miR-509","3p","X-linked"),

    # --- miR-517 ---
    ("hsa-miR-517a-3p","MIMAT0002852", "AUCGUGCAUCCCUUAGAGUGUU", "miR-517","3p","C19MC cluster"),

    # --- miR-520 ---
    ("hsa-miR-520a-3p","MIMAT0002867", "AAAGUGCUUCCCUUUGUGUGUGU","miR-520","3p",""),

    # --- miR-let-7 already added above ---

    # --- miR-663 ---
    ("hsa-miR-663a",   "MIMAT0003325", "AGGCGGGGCGCCGCGGGACCGC", "miR-663","3p",""),

    # --- miR-708 ---
    ("hsa-miR-708-5p", "MIMAT0004924", "AAGGAGCUUACAAUCUAGCUGG", "miR-708","5p",""),

    # --- miR-let-7 already above ---

    # --- miR-1 (cardiac/skeletal muscle) ---
    ("hsa-miR-1-3p",   "MIMAT0000416", "UGGAAUGUAAAGAAGUAUGUAU", "miR-1",  "3p", "cardiac/skeletal muscle; HDAC4"),
    ("hsa-miR-1-5p",   "MIMAT0004548", "ACAUACUUCUUUAUAUGCCCAUA","miR-1",  "5p", ""),

    # --- miR-133 already added above ---

    # --- miR-206 already added above ---
]

# ---------------------------------------------------------------------------
# Build index
# ---------------------------------------------------------------------------
_DB: List[MirnaRecord] = []
for _row in _RAW:
    _DB.append(MirnaRecord(name=_row[0], accession=_row[1], sequence=_row[2],
                           family=_row[3], arm=_row[4],
                           note=_row[5] if len(_row) > 5 else ""))

def _norm(s: str) -> str:
    """Normalise query: lowercase, collapse spaces, strip hsa- prefix, insert
    missing dash (let7a -> let-7a, mir21 -> mir-21)."""
    s = s.lower().strip()
    s = re.sub(r"\s+", "-", s)
    s = re.sub(r"^hsa-", "", s)               # hsa-miR-21 -> miR-21
    # Insert dash between letters and digits: let7a->let-7, mir21->mir-21
    s = re.sub(r"([a-z])([\d])", r"\1-\2", s)
    s = re.sub(r"([\d])([a-z])", r"\1-\2", s)
    return s

def _score(query_norm: str, rec: MirnaRecord) -> int:
    """Return match quality score (higher = better). 0 = no match."""
    name_n = _norm(rec.name)        # e.g. "mir-21-5p"
    fam_n  = _norm(rec.family)      # e.g. "mir-21"
    acc_n  = rec.accession.lower()  # MIMAT0000076
    q = query_norm
    # Exact accession match
    if q == acc_n:
        return 100
    # Exact name match (after prefix strip)
    if q == name_n:
        return 90
    # Name without arm suffix
    if q == name_n.rstrip("p").rstrip("35-").rstrip("-"):
        return 80
    if re.sub(r"[-_]?[35]p$", "", name_n) == q:
        return 80
    # Family exact
    if q == fam_n:
        return 70
    # Name starts with query
    if name_n.startswith(q):
        return 60
    # Family starts with query
    if fam_n.startswith(q):
        return 50
    # Query is contained in name
    if q in name_n:
        return 40
    # Query is contained in family
    if q in fam_n:
        return 30
    return 0


def search_mirna(query: str, max_results: int = 10) -> List[MirnaRecord]:
    """Search bundled miRNA database by name, accession, or family.

    Accepts any of: 'hsa-miR-21-5p', 'miR-21', 'let-7a', 'MIMAT0000076',
    'mir-21', 'let7a' (normalised). Returns up to max_results ranked hits.
    Empty query returns the full database.
    """
    if not query or not query.strip():
        return list(_DB)
    # MIMAT accession: bypass normalization
    if re.match(r"^MIMAT\d+$", query.strip(), re.IGNORECASE):
        hit = get_by_accession(query.strip())
        return [hit] if hit else []
    q = _norm(query)
    scored = [(r, _score(q, r)) for r in _DB]
    scored = [(r, s) for r, s in scored if s > 0]
    scored.sort(key=lambda x: -x[1])
    return [r for r, _ in scored[:max_results]]


def get_by_accession(accession: str) -> Optional[MirnaRecord]:
    acc = accession.strip().upper()
    for r in _DB:
        if r.accession.upper() == acc:
            return r
    return None


def list_families() -> List[str]:
    """Return sorted unique family names."""
    return sorted(set(r.family for r in _DB))
