from __future__ import annotations

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.config import settings
from app.routers import specificity_router, thermo_router, engines_router

app = FastAPI(
    title="FreePrimers Backend",
    description=(
        "Real nearest-neighbour thermodynamics (primer3-py / ViennaRNA) "
        "and BLAST-based specificity checking for the FreePrimers primer "
        "design tool."
    ),
    version="0.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=[o.strip() for o in settings.cors_origins.split(",")],
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(thermo_router.router)
app.include_router(specificity_router.router)
app.include_router(engines_router.router)


@app.get("/health", tags=["meta"])
def health() -> dict[str, str]:
    return {"status": "ok"}
