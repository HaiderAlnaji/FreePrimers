FROM python:3.12-slim

# ncbi-blast+ provides blastn / makeblastdb for the local specificity
# backend. ViennaRNA's Python bindings ship as a prebuilt wheel for
# common platforms, but build-essential is kept in case pip needs to
# build from source on an unsupported architecture.
RUN apt-get update && apt-get install -y --no-install-recommends \
    ncbi-blast+ \
    build-essential \
    curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /srv

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app
COPY scripts ./scripts

ENV PRIMERFORGE_BLAST_DB_DIR=/data/blastdb
RUN mkdir -p /data/blastdb

EXPOSE 8000

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
