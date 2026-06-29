# PrimerForge Backend

Real nearest-neighbour thermodynamics and BLAST-based specificity
checking for the PrimerForge primer design tool. This replaces two
things the client-side (browser-only) version cannot do:

1. **Real thermodynamics.** The frontend's Tm/hairpin/dimer numbers
   are fast heuristics so the tool works with zero setup. This
   backend calls **primer3-py** (the same nearest-neighbour core used
   by Primer3 itself) for Tm, hairpin \u0394G, and heterodimer \u0394G with
   proper structure-finding \u2014 and **ViennaRNA** for RNA secondary
   structure questions (e.g. does a pre-miRNA region fold in a way
   that blocks a primer site).
2. **Specificity.** Whether a primer also binds somewhere it
   shouldn't \u2014 a question that needs a real sequence database and
   BLAST, which cannot run in a browser.

## Quick start

```bash
docker compose build
docker compose run --rm primerforge-api bash scripts/build_databases.sh
docker compose up
```

The API is then at `http://localhost:8000`, with interactive docs at
`http://localhost:8000/docs`.

To confirm everything actually works end to end against your running
instance:

```bash
python3 scripts/smoke_test.py http://localhost:8000
```

This hits every endpoint with known sequences and checks both the
response shape and a few sanity bounds (e.g. that two perfect
complements form a strongly negative-\u0394G duplex). Run it with
`--skip-specificity` before building a database, and with
`--database your_db_name --backend local` once you have one. It's the
fastest way to verify the one thing that couldn't be tested while
building this: that `blastn` is actually installed and reachable on
your machine.

Without Docker, on macOS:

```bash
brew install blast              # provides blastn / makeblastdb
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
PRIMERFORGE_BLAST_DB_DIR=./blastdb bash scripts/build_databases.sh
PRIMERFORGE_BLAST_DB_DIR=./blastdb uvicorn app.main:app --reload
```

Without Docker, on Linux:

```bash
sudo apt install ncbi-blast+    # provides blastn / makeblastdb
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
PRIMERFORGE_BLAST_DB_DIR=./blastdb bash scripts/build_databases.sh
PRIMERFORGE_BLAST_DB_DIR=./blastdb uvicorn app.main:app --reload
```

In both cases, the `PRIMERFORGE_BLAST_DB_DIR` value must match between
the `build_databases.sh` call and the `uvicorn` call — it's just a
folder on disk, name it whatever you like. Run
`source venv/bin/activate` again in any new terminal tab before
running `uvicorn` or the smoke test, since the venv is only active in
the shell session that activated it.

## Endpoints

### Thermodynamics (`/thermo`) \u2014 always available, no setup needed

| Endpoint | Purpose |
|---|---|
| `POST /thermo/tm` | Real NN Tm for one oligo |
| `POST /thermo/hairpin` | Hairpin \u0394G with structure-finding |
| `POST /thermo/duplex` | Heterodimer \u0394G between two oligos (pass the same sequence twice for self-dimer) |
| `POST /thermo/rna-fold` | ViennaRNA MFE structure for an RNA sequence |
| `POST /thermo/panel-screen` | All-pairs cross-reactivity across N named oligos \u2014 use this for an ARMS-PCR panel (checks all 4 primers against each other, not just self-dimers) |

All thermo endpoints accept an optional `salts` object
(`na_mm`, `mg_mm`, `dntp_mm`, `primer_nm`) matching the Advanced panel
in the frontend tool, so results line up with whatever conditions the
user set there.

### Specificity (`POST /specificity`) \u2014 needs a database

```json
{
  "seq": "AGCCTAGGATCCGATTATC",
  "database": "mirbase_mature",
  "backend": "auto",
  "max_hits": 25
}
```

`backend`:
- `"local"` \u2014 blastn against a local database built by
  `scripts/build_databases.sh`. Fast, unlimited, offline. **Use this
  for any production or batch use.**
- `"ncbi"` \u2014 the public NCBI BLAST Common URL API (Put/Get). No
  setup, but each search can take 30s\u2013several minutes and is
  rate-limited by NCBI. Fine for occasional interactive checks only.
- `"auto"` (default) \u2014 uses local if that database name exists,
  otherwise falls back to NCBI and says so in the response's
  `warnings` field.

`GET /specificity/databases` lists local databases currently
available to the server.

Each hit includes percent identity, alignment length, mismatches,
e-value, bit score, and a `risk` classification (`high`/`moderate`/
`low`) based on coverage, identity, and whether the alignment reaches
the primer's 3' end (the end that has to extend for amplification).

## What's verified vs. not

Built and tested in a sandboxed environment with network egress
restricted to a small allowlist. Here's exactly what that means for
trust in this code:

- **`/thermo/*` endpoints: fully tested.** primer3-py and ViennaRNA
  were installed and exercised directly \u2014 Tm, hairpin, heterodimer,
  RNA folding, and the all-pairs panel screen were all run against
  real sequences with sane, expected output (see the duplex example
  above: two perfectly complementary 19-20mers gave \u221218.76 kcal/mol,
  which is what you'd want a real cross-dimer check to catch).
- **`/specificity` local backend: code is correct by construction,
  not live-tested end to end**, because no `blastn` binary exists in
  this sandbox and apt isn't reachable here. The command
  construction, FASTA piping, and tabular-output parsing follow
  standard, well-documented BLAST+ usage \u2014 but you should run
  `scripts/build_databases.sh` and try a real query before relying on
  it for anything important. If something doesn't parse correctly,
  the most likely culprit is a BLAST+ version difference in the
  `-outfmt` column order.
- **`/specificity` NCBI backend: request/response contract verified
  against NCBI's published documentation, not live-tested**, because
  `blast.ncbi.nlm.nih.gov` is outside this sandbox's network
  allowlist (confirmed: requests get a `403 host_not_allowed` from
  the sandbox's own egress proxy, not from NCBI). The Put/Get/RID
  polling flow and parameter names match NCBI's Common URL API
  documentation exactly, but the actual round trip \u2014 in particular
  the regex that extracts hit rows from the Text+Tabular report \u2014
  has not been run against a real NCBI response. **Test this against
  a real query before depending on it**, and treat the local backend
  as the primary, trusted path.
- **`scripts/build_databases.sh`: URL not live-tested** for the same
  network reason. The miRBase download path was corrected after an
  initial wrong guess by checking current third-party pipeline docs;
  if it 404s, check `https://www.mirbase.org/download/` for the
  current location.

## Known simplifications worth knowing about

- The specificity `risk` classification is a heuristic over BLAST's
  summary stats (identity / coverage / whether the 3' end aligns),
  not a true thermodynamic off-target amplification prediction. A
  "low" risk hit could still matter in an unusual reaction; a "high"
  risk hit might not actually amplify under your specific cycling
  conditions. Use it to prioritize what to look at, not as a pass/fail
  gate.
- `panel-screen` checks pairwise duplex formation, not three- or
  four-primer simultaneous interactions, which is a reasonable
  approximation but not exhaustive for complex multiplex panels.
