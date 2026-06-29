# FreePrimers — Backend

Real nearest-neighbour thermodynamics, six assay-specific design engines,
and BLAST-based specificity checking for the FreePrimers primer design platform.

## Quick start (Docker — recommended)

```bash
docker compose build
docker compose up
```

The API is then at `http://localhost:8000`, with interactive docs at
`http://localhost:8000/docs`.

## Quick start (without Docker)

**macOS:**
```bash
brew install blast
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
PYTHONPATH=. uvicorn app.main:app --reload --port 8000
```

**Linux:**
```bash
sudo apt install ncbi-blast+
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
PYTHONPATH=. uvicorn app.main:app --reload --port 8000
```

Run the smoke test to verify everything works:
```bash
python3 scripts/smoke_test.py http://localhost:8000
```

## Endpoints

### Thermodynamics (`/thermo`) — no setup needed

| Endpoint | Purpose |
|---|---|
| `POST /thermo/tm` | Nearest-neighbour Tm for one oligo |
| `POST /thermo/hairpin` | Hairpin ΔG with structure-finding |
| `POST /thermo/duplex` | Heterodimer ΔG between two oligos |
| `POST /thermo/rna-fold` | ViennaRNA MFE structure for an RNA sequence |
| `POST /thermo/panel-screen` | All-pairs cross-reactivity across N oligos |

### Design engines (`/engines`)

| Endpoint | Purpose |
|---|---|
| `POST /engines/arms/design` | Tetra-primer ARMS-PCR optimizer |
| `POST /engines/mirna/design` | Stem-loop miRNA RT-qPCR designer |
| `POST /engines/snp/design` | SNP-aware primer placement |
| `POST /engines/methylation/design` | MSP / BSP / MS-HRM designer |
| `POST /engines/hrm/design` | HRM genotyping (standard + allele-specific) |
| `POST /engines/util/fetch-cdna` | Fetch spliced cDNA from Ensembl |
| `POST /engines/util/fetch-region` | Fetch genomic region from Ensembl |
| `POST /engines/snp/fetch-rsid` | Fetch SNP annotation from Ensembl VEP |

### Specificity (`POST /specificity`) — needs BLAST

```json
{
  "seq": "AGCCTAGGATCCGATTATC",
  "database": "nt",
  "backend": "ncbi",
  "max_hits": 25
}
```

`backend` options:
- `"local"` — blastn against a local database (fast, offline, recommended for batch use)
- `"ncbi"` — NCBI remote BLAST (no setup, but slow and rate-limited)
- `"auto"` — tries local first, falls back to NCBI

## Citation

If you use FreePrimers in your research, please cite:

> Al-Naji H. et al. FreePrimers: an integrated, assay-aware primer design
> platform. *BMC Bioinformatics* (under review, 2025).
> https://github.com/HaiderAlnaji/FreePrimers

## License

[MIT](../LICENSE) © 2025 Haider Al-Naji and contributors
