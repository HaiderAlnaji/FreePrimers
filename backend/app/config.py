"""
Configuration via environment variables (twelve-factor style), so the
same image runs in dev, on your machine, or wherever you deploy it
later without code changes.
"""

from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass
class Settings:
    blast_db_dir: str = os.environ.get("PRIMERFORGE_BLAST_DB_DIR", "/data/blastdb")
    blastn_path: str = os.environ.get("PRIMERFORGE_BLASTN_PATH", "blastn")
    makeblastdb_path: str = os.environ.get("PRIMERFORGE_MAKEBLASTDB_PATH", "makeblastdb")
    cors_origins: str = os.environ.get("PRIMERFORGE_CORS_ORIGINS", "*")


settings = Settings()
