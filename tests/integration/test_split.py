"""Integration tests for spanseq split.

These tests run the full split pipeline using real tools (kma, ccphylo)
and real FASTA data from the repo's test/ directory.

Run with:
    pytest tests/integration/test_split.py -v
"""
import shutil
import pytest
from pathlib import Path

from spanseq.config import SpanSeqConfig
from spanseq.pipeline import SpanSeqPipeline

from tests.integration.conftest import fasta_ids, tsv_column

pytestmark = pytest.mark.integration


# ── Helpers ───────────────────────────────────────────────────────────────────

def _make_cfg(fsa: Path, out: Path, **kwargs) -> SpanSeqConfig:
    """Build a minimal split config for a nucleotide FASTA."""
    defaults = dict(
        action="split",
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

class TestSplitOutputFiles:
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


class TestSplitNoKmerSizeRequired:
    """Regression: kmer_size must be optional (was required in legacy code)."""

    def test_split_works_without_kmer_size(self, macrolide_fsa, tmp_path):
        """Omitting -k should not raise an error."""
        cfg = _make_cfg(macrolide_fsa, tmp_path, kmer_size=None)
        outputs = SpanSeqPipeline(cfg).run()
        assert outputs["makespan"].exists()

    def test_split_works_with_explicit_kmer_size(self, macrolide_fsa, tmp_path):
        """Providing -k 16 should also work."""
        cfg = _make_cfg(macrolide_fsa, tmp_path, kmer_size=16)
        outputs = SpanSeqPipeline(cfg).run()
        assert outputs["makespan"].exists()


class TestSplitSequenceAssignment:
    """Every input sequence must appear in the cluster output."""

    def test_all_sequences_in_clusters(self, macrolide_fsa, tmp_path):
        cfg = _make_cfg(macrolide_fsa, tmp_path)
        outputs = SpanSeqPipeline(cfg).run()

        input_ids = fasta_ids(macrolide_fsa)
        cluster_ids = set(tsv_column(outputs["clusters"], "#Sample"))
        assert input_ids == cluster_ids, (
            f"Sequences missing from clusters: {input_ids - cluster_ids}"
        )

    def test_makespan_covers_all_clusters(self, macrolide_fsa, tmp_path):
        """Every cluster from DBSCAN must appear in the makespan output."""
        cfg = _make_cfg(macrolide_fsa, tmp_path)
        outputs = SpanSeqPipeline(cfg).run()

        cluster_ids = set(tsv_column(outputs["clusters"], "Cluster"))
        makespan_cluster_ids = set(tsv_column(outputs["makespan"], "#Cluster"))
        assert cluster_ids == makespan_cluster_ids, (
            f"Clusters missing from makespan: {cluster_ids - makespan_cluster_ids}"
        )


class TestSplitBins:
    """Verify bin counts and proportional bins."""

    def test_integer_bins(self, macrolide_fsa, tmp_path):
        cfg = _make_cfg(macrolide_fsa, tmp_path, bins=3)
        outputs = SpanSeqPipeline(cfg).run()

        bin_values = set(tsv_column(outputs["makespan"], "Partition"))
        assert len(bin_values) <= 3, f"Expected at most 3 bins, got {bin_values}"

    def test_proportional_bins(self, macrolide_fsa, tmp_path):
        cfg = _make_cfg(macrolide_fsa, tmp_path, bins=[6, 2, 2])
        outputs = SpanSeqPipeline(cfg).run()
        assert outputs["makespan"].exists()

    def test_single_bin(self, macrolide_fsa, tmp_path):
        """With bins=1, all sequences go into one partition."""
        cfg = _make_cfg(macrolide_fsa, tmp_path, bins=1)
        outputs = SpanSeqPipeline(cfg).run()

        bin_values = set(tsv_column(outputs["makespan"], "Partition"))
        assert len(bin_values) == 1


class TestSplitOutputFormats:
    """Test output_format options."""

    def test_merged_table(self, macrolide_fsa, tmp_path):
        cfg = _make_cfg(macrolide_fsa, tmp_path, output_format="merged_table")
        outputs = SpanSeqPipeline(cfg).run()

        assert "partitions" in outputs
        assert outputs["partitions"].exists()

    def test_fasta_files(self, macrolide_fsa, tmp_path):
        cfg = _make_cfg(macrolide_fsa, tmp_path, bins=3, output_format="fasta_files")
        outputs = SpanSeqPipeline(cfg).run()

        assert "fasta_files" in outputs
        assert len(outputs["fasta_files"]) == 3
        for f in outputs["fasta_files"]:
            assert f.exists(), f"Partition FASTA missing: {f}"


class TestSplitDistanceMethods:
    """Smoke tests for each distance method (using the smaller macrolide set)."""

    @pytest.mark.parametrize("method", ["cosine", "jaccard", "szymkiewicz_simpson"])
    def test_kma_methods(self, macrolide_fsa, tmp_path, method):
        cfg = _make_cfg(macrolide_fsa, tmp_path / method, distance_method=method)
        outputs = SpanSeqPipeline(cfg).run()
        assert outputs["makespan"].exists()

    @pytest.mark.skipif(shutil.which("mash") is None, reason="mash not in PATH")
    @pytest.mark.xfail(reason="ccphylo makespan segfaults on mash distance matrices — known tool issue")
    def test_mash(self, macrolide_fsa, tmp_path):
        cfg = _make_cfg(macrolide_fsa, tmp_path, distance_method="mash", max_length=2000)
        outputs = SpanSeqPipeline(cfg).run()
        assert outputs["makespan"].exists()

    @pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs not in PATH")
    def test_mmseqs_fast(self, macrolide_fsa, tmp_path):
        cfg = _make_cfg(macrolide_fsa, tmp_path, distance_method="mmseqs-fast")
        outputs = SpanSeqPipeline(cfg).run()
        assert outputs["makespan"].exists()
