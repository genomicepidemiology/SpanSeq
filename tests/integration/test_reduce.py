"""Integration tests for spanseq reduce.

These tests run the full reduce pipeline using real tools (kma, ccphylo)
and real FASTA data from the repo's test/ directory.

Run with:
    pytest tests/integration/test_reduce.py -v
"""
import pytest
from pathlib import Path

from spanseq.config import SpanSeqConfig
from spanseq.pipeline import SpanSeqPipeline

from tests.integration.conftest import fasta_ids, tsv_column, tsv_column_by_index

pytestmark = pytest.mark.integration


# ── Helpers ───────────────────────────────────────────────────────────────────

def _make_cfg(fsa: Path, out: Path, **kwargs) -> SpanSeqConfig:
    """Build a minimal reduce config for a nucleotide FASTA."""
    defaults = dict(
        action="reduce",
        input_path=fsa,
        input_format="file",
        seq_type="nucleotides",
        output_dir=out,
        min_dist=0.3,
        bins=3,
    )
    defaults.update(kwargs)
    return SpanSeqConfig(**defaults)


# ── Tests ─────────────────────────────────────────────────────────────────────

class TestReduceOutputFiles:
    """Verify that the expected output files are created."""

    def test_creates_clusters_makespan_stats(self, macrolide_fsa, tmp_path):
        cfg = _make_cfg(macrolide_fsa, tmp_path)
        outputs = SpanSeqPipeline(cfg).run()

        assert "clusters" in outputs
        assert "makespan" in outputs
        assert "stats" in outputs
        assert outputs["clusters"].exists(), "clusters TSV missing"
        assert outputs["makespan"].exists(), "makespan TSV missing"
        assert outputs["stats"].exists(),    "stats TSV missing"

    def test_output_files_are_nonempty(self, macrolide_fsa, tmp_path):
        cfg = _make_cfg(macrolide_fsa, tmp_path)
        outputs = SpanSeqPipeline(cfg).run()

        for key in ("clusters", "makespan", "stats"):
            size = outputs[key].stat().st_size
            assert size > 0, f"{key} file is empty"


class TestReduceNoKmerSizeRequired:
    """Regression: kmer_size must be optional (was required in legacy code)."""

    def test_reduce_works_without_kmer_size(self, macrolide_fsa, tmp_path):
        """Omitting -k should not raise an error."""
        cfg = _make_cfg(macrolide_fsa, tmp_path, kmer_size=None)
        outputs = SpanSeqPipeline(cfg).run()
        assert outputs["makespan"].exists()

    def test_reduce_works_with_explicit_kmer_size(self, macrolide_fsa, tmp_path):
        """Providing -k 16 should also work."""
        cfg = _make_cfg(macrolide_fsa, tmp_path, kmer_size=16)
        outputs = SpanSeqPipeline(cfg).run()
        assert outputs["makespan"].exists()


class TestReduceSequenceAssignment:
    """All input sequences must be represented in the cluster output."""

    def test_all_sequences_in_clusters(self, macrolide_fsa, tmp_path):
        cfg = _make_cfg(macrolide_fsa, tmp_path)
        outputs = SpanSeqPipeline(cfg).run()

        input_ids = fasta_ids(macrolide_fsa)
        # KMA hobohm1 output has no header; col 0 = sequence ID
        cluster_ids = set(tsv_column_by_index(outputs["clusters"], 0))
        assert cluster_ids == input_ids, (
            f"Sequences missing from clusters: {input_ids - cluster_ids}"
        )

    def test_makespan_covers_all_clusters(self, macrolide_fsa, tmp_path):
        """Every cluster from the reduce step must appear in the makespan output."""
        cfg = _make_cfg(macrolide_fsa, tmp_path)
        outputs = SpanSeqPipeline(cfg).run()

        # KMA hobohm1 output has no header; col 1 = cluster ID (integer)
        cluster_ids = set(tsv_column_by_index(outputs["clusters"], 1))
        makespan_cluster_ids = set(tsv_column(outputs["makespan"], "#Cluster"))
        assert cluster_ids == makespan_cluster_ids, (
            f"Clusters missing from makespan: {cluster_ids - makespan_cluster_ids}"
        )


class TestReduceBins:
    """Verify bin count behaviour."""

    def test_integer_bins(self, macrolide_fsa, tmp_path):
        cfg = _make_cfg(macrolide_fsa, tmp_path, bins=3)
        outputs = SpanSeqPipeline(cfg).run()

        bin_values = set(tsv_column(outputs["makespan"], "Partition"))
        assert len(bin_values) <= 3, f"Expected at most 3 bins, got {bin_values}"

    def test_single_bin(self, macrolide_fsa, tmp_path):
        cfg = _make_cfg(macrolide_fsa, tmp_path, bins=1)
        outputs = SpanSeqPipeline(cfg).run()

        bin_values = set(tsv_column(outputs["makespan"], "Partition"))
        assert len(bin_values) == 1
