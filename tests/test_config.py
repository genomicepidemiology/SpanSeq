"""Tests for SpanSeqConfig and related utilities."""
import argparse
import pytest
from pathlib import Path

from spanseq.config import (
    SpanSeqConfig,
    compute_sketch_size,
    DISTANCE_TOOLS,
    KMA_DISTANCE_METHODS,
)


class TestComputeSketchSize:
    def test_power_of_two(self):
        assert compute_sketch_size(1024) == 1024

    def test_rounds_up(self):
        assert compute_sketch_size(1000) == 1024

    def test_small(self):
        assert compute_sketch_size(1) == 1
        assert compute_sketch_size(3) == 4

    def test_zero_raises(self):
        with pytest.raises(ValueError):
            compute_sketch_size(0)

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            compute_sketch_size(-5)


class TestDistanceMappings:
    def test_kma_methods(self):
        for method in ("jaccard", "szymkiewicz_simpson", "cosine", "kmer_inv"):
            assert DISTANCE_TOOLS[method] == "kma"

    def test_mash(self):
        assert DISTANCE_TOOLS["mash"] == "mash"

    def test_identity(self):
        assert DISTANCE_TOOLS["identity"] == "ggsearch36"

    def test_mmseqs2(self):
        assert DISTANCE_TOOLS["mmseqs2"] == "mmseqs2"


class TestSpanSeqConfig:
    def _make_config(self, tmp_path, **overrides):
        defaults = dict(
            action="split",
            input_path=tmp_path / "input.fsa",
            input_format="file",
            seq_type="nucleotides",
            output_dir=tmp_path / "output",
        )
        # Create the input file so Path.resolve() works
        defaults["input_path"].parent.mkdir(parents=True, exist_ok=True)
        defaults["input_path"].write_text(">s\nA\n")
        defaults.update(overrides)
        return SpanSeqConfig(**defaults)

    def test_paths_resolved(self, tmp_path):
        cfg = self._make_config(tmp_path)
        assert cfg.input_path.is_absolute()
        assert cfg.output_dir.is_absolute()

    def test_tmp_dir_default(self, tmp_path):
        cfg = self._make_config(tmp_path)
        assert cfg.tmp_dir == cfg.output_dir / "tmp"

    def test_tmp_dir_explicit(self, tmp_path):
        cfg = self._make_config(tmp_path, tmp_dir=tmp_path / "custom_tmp")
        assert cfg.tmp_dir == (tmp_path / "custom_tmp").resolve()

    def test_distance_tool_cosine(self, tmp_path):
        cfg = self._make_config(tmp_path, distance_method="cosine")
        assert cfg.distance_tool == "kma"

    def test_distance_tool_mash(self, tmp_path):
        cfg = self._make_config(tmp_path, distance_method="mash")
        assert cfg.distance_tool == "mash"

    def test_distance_tool_identity(self, tmp_path):
        cfg = self._make_config(tmp_path, distance_method="identity")
        assert cfg.distance_tool == "ggsearch36"

    def test_distance_tool_mmseqs2(self, tmp_path):
        cfg = self._make_config(tmp_path, distance_method="mmseqs2")
        assert cfg.distance_tool == "mmseqs2"

    def test_kma_dist_flag(self, tmp_path):
        cfg = self._make_config(tmp_path, distance_method="cosine")
        assert cfg.kma_dist_flag == 256

    def test_kma_dist_flag_jaccard(self, tmp_path):
        cfg = self._make_config(tmp_path, distance_method="jaccard")
        assert cfg.kma_dist_flag == 64

    def test_effective_dist_value(self, tmp_path):
        cfg = self._make_config(tmp_path, distance_method="cosine", min_dist=0.5)
        assert cfg.effective_dist_value == 0.5  # factor is 1.0

    def test_effective_dist_kmer_inv(self, tmp_path):
        cfg = self._make_config(tmp_path, distance_method="kmer_inv", min_dist=50.0)
        assert cfg.effective_dist_value == pytest.approx(0.5)  # factor is 100.0

    def test_needs_hobohm(self, tmp_path):
        cfg = self._make_config(tmp_path, approach="hobohm_reduce")
        assert cfg.needs_hobohm is True

    def test_not_needs_hobohm(self, tmp_path):
        cfg = self._make_config(tmp_path, approach="all")
        assert cfg.needs_hobohm is False

    def test_sketch_size_computed_for_mash(self, tmp_path):
        cfg = self._make_config(tmp_path, distance_method="mash", max_length=1000)
        assert cfg.sketch_size == 1024

    def test_sample_name(self, tmp_path):
        cfg = self._make_config(tmp_path)
        assert cfg.sample_name == "input"


class TestFromArgs:
    def _make_args(self, **overrides):
        defaults = dict(
            action="split",
            input_fasta="/tmp/test.fsa",
            input_folder=None,
            input_batch=None,
            seqtype="nucleotides",
            output_folder="/tmp/out",
            output_format="minimal",
            keep_tmp=False,
            min_dist=0.3,
            bins="5",
            makespanProcess="DBF",
            makespanWeights="none",
            kmer_size=None,
            minimizer_size=None,
            prefix="-",
            MegaDB=False,
            max_length=None,
            threads=1,
            memory_disk=False,
            temp_files=None,
            distanceMethod="cosine",
            approach="all",
            hobohm1_distance=None,
            hobohm1_method="cdhit",
            makespanImbalanced=None,
            kmaPath=None,
            mashPath=None,
            CDHitPath=None,
            ccphyloPath=None,
            GGSearchPath=None,
            mmseqsPath=None,
        )
        defaults.update(overrides)
        return argparse.Namespace(**defaults)

    def test_from_args_basic(self):
        args = self._make_args()
        cfg = SpanSeqConfig.from_args(args)
        assert cfg.action == "split"
        assert cfg.input_format == "file"
        assert cfg.distance_method == "cosine"

    def test_from_args_reduce(self):
        args = self._make_args(action="reduce")
        cfg = SpanSeqConfig.from_args(args)
        assert cfg.action == "reduce"

    def test_from_args_folder_input(self):
        args = self._make_args(input_fasta=None, input_folder="/tmp/seqs")
        cfg = SpanSeqConfig.from_args(args)
        assert cfg.input_format == "folder"

    def test_from_args_bins_csv(self):
        args = self._make_args(bins="3,4,5")
        cfg = SpanSeqConfig.from_args(args)
        assert cfg.bins == [3, 4, 5]

    def test_from_args_bins_int(self):
        args = self._make_args(bins="5")
        cfg = SpanSeqConfig.from_args(args)
        assert cfg.bins == 5

    def test_from_args_no_input_raises(self):
        args = self._make_args(input_fasta=None, input_folder=None, input_batch=None)
        with pytest.raises(ValueError, match="No input"):
            SpanSeqConfig.from_args(args)
