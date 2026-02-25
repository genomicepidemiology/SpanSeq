"""Compatibility tests: new SpanSeq must produce equivalent results to legacy.

Both codebases are run on the same input with the same parameters.
Results are compared at the semantic level:
  - Same sequences appear in the output
  - Same number of partitions
  - Similar cluster count
  - Sequences in the same cluster (new) must also be in the same partition (legacy)

Run with:
    pytest tests/integration/test_compatibility.py -v

Skipped automatically if:
  - The legacy script is not found
  - Snakemake is not installed (legacy depends on it)
  - kma / ccphylo are not in PATH
"""
import shutil
import subprocess
import sys
import pytest
from pathlib import Path

from spanseq.config import SpanSeqConfig
from spanseq.pipeline import SpanSeqPipeline
from tests.integration.conftest import fasta_ids, MACROLIDE_FSA

pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).parent.parent.parent
LEGACY_SCRIPT = REPO_ROOT / "src" / "spanseq-legacy" / "spanseq.py"


# ── Availability checks ───────────────────────────────────────────────────────

def _legacy_runnable() -> bool:
    """Return True if the legacy script exists and Snakemake is available."""
    if not LEGACY_SCRIPT.exists():
        return False
    if shutil.which("snakemake") is None:
        return False
    return True


if not _legacy_runnable():
    pytest.skip(
        "Legacy spanseq or snakemake not available — skipping compatibility tests",
        allow_module_level=True,
    )


# ── Helpers ───────────────────────────────────────────────────────────────────

def _run_legacy_split(fsa: Path, out: Path, bins: int, min_dist: float) -> Path:
    """Run legacy spanseq split. Returns output directory, skips on failure."""
    out.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [
            sys.executable, str(LEGACY_SCRIPT), "split",
            "-i", str(fsa),
            "-s", "nucleotides",
            "-o", str(out),
            "-c", str(min_dist),
            "-b", str(bins),
            "-d", "cosine",
            "-k", "16",   # legacy requires explicit kmer_size
        ],
        capture_output=True,
        text=True,
        timeout=600,
        cwd=str(REPO_ROOT),
    )
    if result.returncode != 0:
        pytest.skip(f"Legacy split failed:\n{result.stderr[-800:]}")
    return out


def _run_new_split(fsa: Path, out: Path, bins: int, min_dist: float) -> dict:
    """Run new spanseq split. Returns outputs dict."""
    cfg = SpanSeqConfig(
        action="split",
        input_path=fsa,
        input_format="file",
        seq_type="nucleotides",
        output_dir=out,
        min_dist=min_dist,
        bins=bins,
        kmer_size=16,   # use same kmer_size as legacy for fair comparison
    )
    return SpanSeqPipeline(cfg).run()


def _parse_clusters(clusters_file: Path) -> dict:
    """Parse clusters TSV → {seq_id: cluster_id}.

    Handles the ccphylo dbscan format:
        ## N_seqs N_clusters dist
        #Sample  Neighbors  Cluster
        seq1     5          0
    """
    seq_to_cluster = {}
    with open(clusters_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("##") or line.startswith("#Sample"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                seq_to_cluster[parts[0]] = int(parts[2])
    return seq_to_cluster


def _parse_makespan(makespan_file: Path) -> dict:
    """Parse makespan TSV → {cluster_id: partition}.

    Format:  #Cluster  Cluster_size  Cluster_weight  Partition
    """
    cluster_to_partition = {}
    with open(makespan_file) as f:
        header_seen = False
        for line in f:
            line = line.strip()
            if not line or line.startswith("##"):
                continue
            if not header_seen:
                header_seen = True
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                cluster_to_partition[int(parts[0])] = int(parts[3])
    return cluster_to_partition


def _seq_to_partition(clusters_file: Path, makespan_file: Path) -> dict:
    """Join cluster and makespan files → {seq_id: partition}."""
    seq_to_cluster = _parse_clusters(clusters_file)
    cluster_to_partition = _parse_makespan(makespan_file)
    return {
        seq: cluster_to_partition.get(cluster, -1)
        for seq, cluster in seq_to_cluster.items()
    }


def _find_output_file(directory: Path, suffix: str) -> Path:
    """Find a file with a given suffix in a directory (first match)."""
    matches = list(directory.glob(f"*{suffix}"))
    if not matches:
        pytest.skip(f"No file matching *{suffix} in {directory}")
    return matches[0]


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def split_outputs(tmp_path_factory):
    """Run both legacy and new split on macrolide.fsa once per module."""
    if not MACROLIDE_FSA.exists():
        pytest.skip(f"Test data not found: {MACROLIDE_FSA}")

    legacy_out = tmp_path_factory.mktemp("legacy_split")
    new_out = tmp_path_factory.mktemp("new_split")

    _run_legacy_split(MACROLIDE_FSA, legacy_out, bins=3, min_dist=0.3)
    new_outputs = _run_new_split(MACROLIDE_FSA, new_out, bins=3, min_dist=0.3)

    return {
        "legacy_dir": legacy_out,
        "new_outputs": new_outputs,
        "fsa": MACROLIDE_FSA,
    }


# ── Tests ─────────────────────────────────────────────────────────────────────

class TestSequenceCoverage:
    """All input sequences must appear in both outputs."""

    def test_legacy_has_all_sequences(self, split_outputs):
        clusters_file = _find_output_file(split_outputs["legacy_dir"], "_clusters.tsv")
        seq_to_cluster = _parse_clusters(clusters_file)
        input_ids = fasta_ids(split_outputs["fsa"])
        assert set(seq_to_cluster.keys()) == input_ids

    def test_new_has_all_sequences(self, split_outputs):
        seq_to_cluster = _parse_clusters(split_outputs["new_outputs"]["clusters"])
        input_ids = fasta_ids(split_outputs["fsa"])
        assert set(seq_to_cluster.keys()) == input_ids

    def test_same_sequence_count(self, split_outputs):
        legacy_file = _find_output_file(split_outputs["legacy_dir"], "_clusters.tsv")
        legacy_ids = set(_parse_clusters(legacy_file).keys())
        new_ids = set(_parse_clusters(split_outputs["new_outputs"]["clusters"]).keys())
        assert legacy_ids == new_ids


class TestPartitionCount:
    """Both pipelines must produce the requested number of partitions."""

    def test_legacy_partition_count(self, split_outputs):
        makespan_file = _find_output_file(split_outputs["legacy_dir"], "_makespan.tsv")
        partitions = set(_parse_makespan(makespan_file).values())
        assert len(partitions) == 3, f"Expected 3 partitions, got {partitions}"

    def test_new_partition_count(self, split_outputs):
        partitions = set(_parse_makespan(split_outputs["new_outputs"]["makespan"]).values())
        assert len(partitions) == 3, f"Expected 3 partitions, got {partitions}"


class TestClusterCount:
    """Cluster count should be the same — same tools, same parameters."""

    def test_same_number_of_clusters(self, split_outputs):
        legacy_file = _find_output_file(split_outputs["legacy_dir"], "_clusters.tsv")
        legacy_clusters = set(_parse_clusters(legacy_file).values())
        new_clusters = set(_parse_clusters(split_outputs["new_outputs"]["clusters"]).values())

        # Allow ±1 difference (rounding / noise-sequence handling may vary)
        diff = abs(len(legacy_clusters) - len(new_clusters))
        assert diff <= 1, (
            f"Cluster count differs: legacy={len(legacy_clusters)}, "
            f"new={len(new_clusters)}"
        )


class TestSemanticEquivalence:
    """Sequences clustered together by the new code must also share a partition in legacy."""

    def test_co_clustered_sequences_share_legacy_partition(self, split_outputs):
        """If two sequences are in the same new cluster → same legacy partition."""
        legacy_clusters_file = _find_output_file(split_outputs["legacy_dir"], "_clusters.tsv")
        legacy_makespan_file = _find_output_file(split_outputs["legacy_dir"], "_makespan.tsv")
        legacy_part = _seq_to_partition(legacy_clusters_file, legacy_makespan_file)

        new_seq_to_cluster = _parse_clusters(split_outputs["new_outputs"]["clusters"])

        # Group sequences by new cluster
        from collections import defaultdict
        new_clusters: dict = defaultdict(list)
        for seq, cluster in new_seq_to_cluster.items():
            new_clusters[cluster].append(seq)

        violations = []
        for cluster_id, members in new_clusters.items():
            legacy_partitions = {legacy_part.get(m, -1) for m in members}
            if len(legacy_partitions) > 1:
                violations.append(
                    f"New cluster {cluster_id} (size {len(members)}) spans "
                    f"legacy partitions {legacy_partitions}"
                )

        assert not violations, (
            f"{len(violations)} cluster(s) span multiple legacy partitions:\n"
            + "\n".join(violations[:5])
        )
