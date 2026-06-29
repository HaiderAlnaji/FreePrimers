from __future__ import annotations

from fastapi import APIRouter, HTTPException

from app.models.schemas import SpecificityRequest, SpecificityResult
from app.services import specificity

router = APIRouter(prefix="/specificity", tags=["specificity"])


@router.post("", response_model=SpecificityResult)
async def check(req: SpecificityRequest) -> SpecificityResult:
    """
    BLAST a primer/probe sequence against a reference set and report
    where else it binds. backend='auto' (default) uses a local
    database if one is configured under that name, otherwise falls
    back to the public NCBI BLAST API (slow, rate-limited — fine for
    occasional interactive checks, not for batch use).
    """
    try:
        return await specificity.check_specificity(req)
    except specificity.SpecificityError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@router.get("/databases")
def databases() -> dict[str, list[str]]:
    """List local BLAST databases currently available to this server."""
    return {"local_databases": specificity.list_local_databases()}
