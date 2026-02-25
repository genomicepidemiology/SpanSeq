"""Shared fixtures for SpanSeq integration tests.

Integration tests call real external tools (kma, ccphylo, etc.) and use
the FASTA files shipped in the repo's test/ directory.

Run integration tests:
    pytest tests/integration/ -v

Run only unit tests (default):
    pytest
"""
import shutil
import pytest
from pathlib import Path


REPO_ROOT = Path(__file__).parent.parent.parent
TEST_DATA = REPO_ROOT / "test"

MACROLIDE_FSA = TEST_DATA / "macrolide.fsa"      # 174 sequences — fast
BETA_LACTAM_FSA = TEST_DATA / "beta-lactam.fsa"  # 2021 sequences — slower


# ── Tool checks ───────────────────────────────────────────────────────────────

def _missing(*tools: str) -> list:
    return [t for t in tools if shutil.which(t) is None]


@pytest.fixture(scope="session", autouse=True)
def require_kma_ccphylo():
    """Skip the entire suite if kma or ccphylo are not available."""
    missing = _missing("kma", "ccphylo")
    if missing:
        pytest.skip(f"Required tools not in PATH: {', '.join(missing)}")


# ── Data fixtures ─────────────────────────────────────────────────────────────

@pytest.fixture(scope="session")
def macrolide_fsa() -> Path:
    if not MACROLIDE_FSA.exists():
        pytest.skip(f"Test data not found: {MACROLIDE_FSA}")
    return MACROLIDE_FSA


@pytest.fixture(scope="session")
def beta_lactam_fsa() -> Path:
    if not BETA_LACTAM_FSA.exists():
        pytest.skip(f"Test data not found: {BETA_LACTAM_FSA}")
    return BETA_LACTAM_FSA


# ── Helpers ───────────────────────────────────────────────────────────────────

def fasta_ids(path: Path) -> set:
    """Return the set of sequence IDs from a FASTA file."""
    ids = set()
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def tsv_column(path: Path, header: str) -> list:
    """Return all values in a named TSV column.

    Skips ccphylo metadata lines that start with '##' before reading the header.
    """
    with open(path) as f:
        lines = [l for l in f.readlines() if not l.startswith("##")]
    cols = lines[0].strip().split("\t")
    idx = cols.index(header)
    return [line.strip().split("\t")[idx] for line in lines[1:] if line.strip()]


def tsv_column_by_index(path: Path, col_idx: int) -> list:
    """Return all values in a TSV column by position (0-based).

    For headerless files such as KMA hobohm1 cluster output.
    """
    with open(path) as f:
        lines = [l for l in f.readlines() if not l.startswith("##")]
    return [line.strip().split("\t")[col_idx] for line in lines if line.strip()]
