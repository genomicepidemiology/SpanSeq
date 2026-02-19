"""Tests for SpanSeqPipeline orchestration logic."""
import pytest
from pathlib import Path
from unittest.mock import MagicMock, call

from spanseq.config import SpanSeqConfig
from spanseq.fasta import parse_fasta, write_fasta
from spanseq.pipeline import SpanSeqPipeline


def _make_config(tmp_path, **overrides):
    """Create a SpanSeqConfig with sensible defaults for testing."""
    fasta = tmp_path / "input.fsa"
    fasta.write_text(">s1\nATCG\n>s2\nGGCC\n")
    defaults = dict(
        action="split",
        input_path=fasta,
        input_format="file",
        seq_type="nucleotides",
        output_dir=tmp_path / "output",
        distance_method="cosine",
        min_dist=0.3,
        bins=3,
    )
    defaults.update(overrides)
    return SpanSeqConfig(**defaults)


def _make_pipeline(tmp_path, mocker, **overrides):
    """Helper to create a pipeline with mocked tool resolution."""
    fake_exe = tmp_path / "tool"
    fake_exe.write_text("#!/bin/sh\n")
    fake_exe.chmod(0o755)
    mocker.patch.object(
        SpanSeqPipeline, "_resolve_tool", return_value=fake_exe
    )
    cfg = _make_config(tmp_path, **overrides)
    return SpanSeqPipeline(cfg)


# ── Tool resolution ──────────────────────────────────────────────────


class TestResolveToolAndInit:
    def test_resolve_tool_from_which(self, tmp_path, mocker):
        """_resolve_tool falls back to shutil.which."""
        fake_exe = tmp_path / "ccphylo"
        fake_exe.write_text("#!/bin/sh\n")
        fake_exe.chmod(0o755)
        mocker.patch("shutil.which", return_value=str(fake_exe))

        cfg = _make_config(tmp_path)
        mocker.patch.object(
            SpanSeqPipeline, "_resolve_tool", return_value=fake_exe
        )
        pipeline = SpanSeqPipeline(cfg)
        assert pipeline.dbscan is not None
        assert pipeline.makespan is not None

    def test_resolve_tool_not_found(self, tmp_path, mocker):
        """_resolve_tool raises when tool not in PATH."""
        mocker.patch("shutil.which", return_value=None)
        cfg = _make_config(tmp_path)
        with pytest.raises(FileNotFoundError, match="not found in PATH"):
            SpanSeqPipeline(cfg)

    def test_resolve_tool_user_path(self, tmp_path):
        """_resolve_tool uses user-provided path."""
        exe = tmp_path / "my_tool"
        exe.write_text("#!/bin/sh\n")
        exe.chmod(0o755)

        cfg = _make_config(tmp_path)
        pipeline_cls = SpanSeqPipeline.__new__(SpanSeqPipeline)
        pipeline_cls.config = cfg
        result = pipeline_cls._resolve_tool(exe, "my_tool")
        assert result == exe

    def test_resolve_tool_directory_with_candidate(self, tmp_path):
        """_resolve_tool finds executable inside a user-provided directory."""
        tool_dir = tmp_path / "tools"
        tool_dir.mkdir()
        exe = tool_dir / "ccphylo"
        exe.write_text("#!/bin/sh\n")
        exe.chmod(0o755)

        cfg = _make_config(tmp_path)
        pipeline_cls = SpanSeqPipeline.__new__(SpanSeqPipeline)
        pipeline_cls.config = cfg
        result = pipeline_cls._resolve_tool(tool_dir, "ccphylo")
        assert result == exe

    def test_resolve_tool_directory_without_candidate(self, tmp_path):
        """_resolve_tool returns dir/default even if it doesn't exist yet."""
        tool_dir = tmp_path / "tools"
        tool_dir.mkdir()

        cfg = _make_config(tmp_path)
        pipeline_cls = SpanSeqPipeline.__new__(SpanSeqPipeline)
        pipeline_cls.config = cfg
        result = pipeline_cls._resolve_tool(tool_dir, "ccphylo")
        assert result == tool_dir / "ccphylo"

    def test_resolve_tool_via_which(self, tmp_path, mocker):
        """_resolve_tool calls shutil.which when user_path is None."""
        mocker.patch("shutil.which", return_value="/usr/bin/ccphylo")
        cfg = _make_config(tmp_path)
        pipeline_cls = SpanSeqPipeline.__new__(SpanSeqPipeline)
        pipeline_cls.config = cfg
        result = pipeline_cls._resolve_tool(None, "ccphylo")
        assert result == Path("/usr/bin/ccphylo")


# ── Init apps for different configs ──────────────────────────────────


class TestInitApps:
    def test_init_reduce(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, action="reduce")
        assert hasattr(pipeline, "kma_index")

    def test_init_split_mash(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="mash")
        assert hasattr(pipeline, "mash")

    def test_init_split_ggsearch(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="identity")
        assert hasattr(pipeline, "ggsearch")

    def test_init_split_mmseqs2(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="mmseqs2")
        assert hasattr(pipeline, "mmseqs2")

    def test_init_split_mmseqs_fast(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="mmseqs-fast")
        assert hasattr(pipeline, "mmseqs2_cluster")

    def test_init_split_tree(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, tree=True)
        assert hasattr(pipeline, "tree_app")

    def test_init_split_hobohm_cdhit(self, tmp_path, mocker):
        pipeline = _make_pipeline(
            tmp_path, mocker,
            approach="hobohm_reduce",
            hobohm1_distance=0.05,
            hobohm1_method="cdhit",
        )
        assert hasattr(pipeline, "cdhit")

    def test_init_split_hobohm_kma(self, tmp_path, mocker):
        pipeline = _make_pipeline(
            tmp_path, mocker,
            approach="hobohm_split",
            hobohm1_distance=0.05,
            hobohm1_method="kma",
        )
        assert hasattr(pipeline, "kma_index")


# ── Distance dispatch ────────────────────────────────────────────────


class TestComputeDistanceDispatch:
    def test_dispatch_kma(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="cosine")
        mock_kma = mocker.patch.object(pipeline, "_distance_kma", return_value=Path("x"))
        pipeline._compute_distance(Path("in.fsa"), Path("out.phy"))
        mock_kma.assert_called_once()

    def test_dispatch_mash(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="mash")
        mock_mash = mocker.patch.object(pipeline, "_distance_mash", return_value=Path("x"))
        pipeline._compute_distance(Path("in.fsa"), Path("out.phy"))
        mock_mash.assert_called_once()

    def test_dispatch_ggsearch(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="identity")
        mock_gg = mocker.patch.object(pipeline, "_distance_ggsearch", return_value=Path("x"))
        pipeline._compute_distance(Path("in.fsa"), Path("out.phy"))
        mock_gg.assert_called_once()

    def test_dispatch_mmseqs2(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="mmseqs2")
        mock_mm = mocker.patch.object(pipeline, "_distance_mmseqs2", return_value=Path("x"))
        pipeline._compute_distance(Path("in.fsa"), Path("out.phy"))
        mock_mm.assert_called_once()


# ── Run dispatch ─────────────────────────────────────────────────────


class TestRunDispatch:
    def test_unknown_action_raises(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker)
        pipeline.config.action = "unknown"
        with pytest.raises(ValueError, match="Unknown action"):
            pipeline.run()

    def test_run_dispatches_split(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker)
        mock_split = mocker.patch.object(pipeline, "_run_split", return_value={})
        pipeline.run()
        mock_split.assert_called_once()

    def test_run_dispatches_reduce(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, action="reduce")
        mock_reduce = mocker.patch.object(pipeline, "_run_reduce", return_value={})
        pipeline.run()
        mock_reduce.assert_called_once()


# ── _run_split pipeline ──────────────────────────────────────────────


class TestRunSplit:
    def test_split_minimal(self, tmp_path, mocker):
        """Minimal split: distance → dbscan → makespan, no output formatting."""
        pipeline = _make_pipeline(tmp_path, mocker)
        cfg = pipeline.config

        # Mock distance computation
        dist_file = cfg.tmp_dir / f"{cfg.sample_name}.phy"
        mocker.patch.object(
            pipeline, "_compute_distance",
            side_effect=lambda input_file, output_file: (output_file.parent.mkdir(parents=True, exist_ok=True), output_file.write_text("matrix"), output_file)[-1],
        )

        # Mock dbscan — write cluster output
        def fake_dbscan_run(cmd, workdir, **kw):
            clusters = cfg.results_dir / f"{cfg.sample_name}_clusters.tsv"
            clusters.write_text("#Sample\tCluster\tWeight\ns1\t1\t1\ns2\t1\t1\n")
        mocker.patch.object(pipeline.dbscan, "run", side_effect=fake_dbscan_run)

        # Mock makespan
        def fake_makespan(**kw):
            out = Path(kw["output_file"])
            out.write_text("#Sample\tPartition\ns1\t1\ns2\t2\n")
            stats = Path(kw["stats_file"])
            stats.write_text("stats\n")
            return stats
        mocker.patch.object(pipeline.makespan, "run_makespan", side_effect=fake_makespan)

        outputs = pipeline.run()
        assert "clusters" in outputs
        assert "makespan" in outputs
        assert "stats" in outputs

    def test_split_merged_table(self, tmp_path, mocker):
        """Split with output_format=merged_table calls _merge_tables."""
        pipeline = _make_pipeline(tmp_path, mocker, output_format="merged_table")
        cfg = pipeline.config

        mocker.patch.object(
            pipeline, "_compute_distance",
            side_effect=lambda input_file, output_file: (output_file.parent.mkdir(parents=True, exist_ok=True), output_file.write_text("m"), output_file)[-1],
        )

        def fake_dbscan_run(cmd, workdir, **kw):
            clusters = cfg.results_dir / f"{cfg.sample_name}_clusters.tsv"
            clusters.write_text("#Sample\tCluster\ns1\t1\ns2\t1\n")
        mocker.patch.object(pipeline.dbscan, "run", side_effect=fake_dbscan_run)

        def fake_makespan(**kw):
            out = Path(kw["output_file"])
            out.write_text("#Sample\tPartition\ns1\t1\ns2\t2\n")
            return Path(kw["stats_file"])
        mocker.patch.object(pipeline.makespan, "run_makespan", side_effect=fake_makespan)

        outputs = pipeline.run()
        assert "partitions" in outputs

    def test_split_fasta_files(self, tmp_path, mocker):
        """Split with output_format=fasta_files calls _create_fasta_partitions."""
        pipeline = _make_pipeline(tmp_path, mocker, output_format="fasta_files", bins=2)
        cfg = pipeline.config

        mocker.patch.object(
            pipeline, "_compute_distance",
            side_effect=lambda input_file, output_file: (output_file.parent.mkdir(parents=True, exist_ok=True), output_file.write_text("m"), output_file)[-1],
        )

        def fake_dbscan_run(cmd, workdir, **kw):
            clusters = cfg.results_dir / f"{cfg.sample_name}_clusters.tsv"
            clusters.write_text("#Sample\tCluster\ns1\t1\ns2\t1\n")
        mocker.patch.object(pipeline.dbscan, "run", side_effect=fake_dbscan_run)

        def fake_makespan(**kw):
            out = Path(kw["output_file"])
            out.write_text("#Sample\tPartition\ns1\t1\ns2\t2\n")
            return Path(kw["stats_file"])
        mocker.patch.object(pipeline.makespan, "run_makespan", side_effect=fake_makespan)

        mock_merge = mocker.patch.object(
            pipeline, "_merge_tables",
            side_effect=lambda clusters_file, makespan_file, output_file: (
                output_file.write_text("id\tpartition\ns1\t1\ns2\t2\n"), output_file
            )[-1],
        )
        mock_fasta = mocker.patch.object(
            pipeline, "_create_fasta_partitions",
            return_value=[Path("f1.fsa"), Path("f2.fsa")],
        )

        outputs = pipeline.run()
        assert "fasta_files" in outputs
        mock_merge.assert_called_once()
        mock_fasta.assert_called_once()

    def test_split_with_hobohm(self, tmp_path, mocker):
        """Split with hobohm_reduce calls _run_hobohm1 before distance."""
        pipeline = _make_pipeline(
            tmp_path, mocker,
            approach="hobohm_reduce",
            hobohm1_distance=0.05,
            hobohm1_method="cdhit",
        )
        cfg = pipeline.config

        hobohm_output = cfg.tmp_dir / "hobohm.fsa"
        mocker.patch.object(
            pipeline, "_run_hobohm1", return_value=hobohm_output,
        )

        mocker.patch.object(
            pipeline, "_compute_distance",
            side_effect=lambda input_file, output_file: (output_file.parent.mkdir(parents=True, exist_ok=True), output_file.write_text("m"), output_file)[-1],
        )

        def fake_dbscan_run(cmd, workdir, **kw):
            clusters = cfg.results_dir / f"{cfg.sample_name}_clusters.tsv"
            clusters.write_text("#Sample\tCluster\ns1\t1\n")
        mocker.patch.object(pipeline.dbscan, "run", side_effect=fake_dbscan_run)

        def fake_makespan(**kw):
            out = Path(kw["output_file"])
            out.write_text("#Sample\tPartition\ns1\t1\n")
            return Path(kw["stats_file"])
        mocker.patch.object(pipeline.makespan, "run_makespan", side_effect=fake_makespan)

        outputs = pipeline.run()
        pipeline._run_hobohm1.assert_called_once()

    def test_split_with_imbalance(self, tmp_path, mocker):
        """Split with imbalance_file calls _add_class_columns."""
        imb = tmp_path / "imb.tsv"
        imb.write_text("s1\tA\ns2\tB\n")
        pipeline = _make_pipeline(tmp_path, mocker, imbalance_file=imb)
        cfg = pipeline.config

        mocker.patch.object(
            pipeline, "_compute_distance",
            side_effect=lambda input_file, output_file: (output_file.parent.mkdir(parents=True, exist_ok=True), output_file.write_text("m"), output_file)[-1],
        )

        def fake_dbscan_run(cmd, workdir, **kw):
            clusters = cfg.results_dir / f"{cfg.sample_name}_clusters.tsv"
            clusters.write_text("#Sample\tCluster\ns1\t1\ns2\t1\n")
        mocker.patch.object(pipeline.dbscan, "run", side_effect=fake_dbscan_run)

        mock_add_class = mocker.patch.object(
            pipeline, "_add_class_columns",
            side_effect=lambda clusters_file, output_file: (output_file.write_text("class_data\n"), output_file)[-1],
        )

        def fake_makespan(**kw):
            out = Path(kw["output_file"])
            out.write_text("#Sample\tPartition\ns1\t1\ns2\t2\n")
            return Path(kw["stats_file"])
        mocker.patch.object(pipeline.makespan, "run_makespan", side_effect=fake_makespan)

        outputs = pipeline.run()
        mock_add_class.assert_called_once()

    def test_split_cleanup(self, tmp_path, mocker):
        """Split with keep_tmp=False removes tmp directory."""
        pipeline = _make_pipeline(tmp_path, mocker, keep_tmp=False)
        cfg = pipeline.config

        mocker.patch.object(
            pipeline, "_compute_distance",
            side_effect=lambda input_file, output_file: (output_file.parent.mkdir(parents=True, exist_ok=True), output_file.write_text("m"), output_file)[-1],
        )

        def fake_dbscan_run(cmd, workdir, **kw):
            clusters = cfg.results_dir / f"{cfg.sample_name}_clusters.tsv"
            clusters.write_text("#Sample\tCluster\ns1\t1\n")
        mocker.patch.object(pipeline.dbscan, "run", side_effect=fake_dbscan_run)

        def fake_makespan(**kw):
            out = Path(kw["output_file"])
            out.write_text("#Sample\tPartition\ns1\t1\n")
            return Path(kw["stats_file"])
        mocker.patch.object(pipeline.makespan, "run_makespan", side_effect=fake_makespan)

        mock_rmtree = mocker.patch("shutil.rmtree")
        pipeline.run()
        mock_rmtree.assert_called_once()

    def test_split_with_tree(self, tmp_path, mocker):
        """Split with --tree builds a Newick tree from the distance matrix."""
        pipeline = _make_pipeline(tmp_path, mocker, tree=True, tree_method="nj")
        cfg = pipeline.config

        mocker.patch.object(
            pipeline, "_compute_distance",
            side_effect=lambda input_file, output_file: (output_file.parent.mkdir(parents=True, exist_ok=True), output_file.write_text("matrix"), output_file)[-1],
        )

        def fake_dbscan_run(cmd, workdir, **kw):
            clusters = cfg.results_dir / f"{cfg.sample_name}_clusters.tsv"
            clusters.write_text("#Sample\tCluster\ns1\t1\ns2\t1\n")
        mocker.patch.object(pipeline.dbscan, "run", side_effect=fake_dbscan_run)

        mock_tree_run = mocker.patch.object(pipeline.tree_app, "run")

        def fake_makespan(**kw):
            out = Path(kw["output_file"])
            out.write_text("#Sample\tPartition\ns1\t1\ns2\t2\n")
            return Path(kw["stats_file"])
        mocker.patch.object(pipeline.makespan, "run_makespan", side_effect=fake_makespan)

        outputs = pipeline.run()
        mock_tree_run.assert_called_once()
        assert "tree" in outputs
        assert str(outputs["tree"]).endswith(".nwk")

    def test_split_mmseqs_fast(self, tmp_path, mocker):
        """mmseqs-fast skips distance+dbscan, calls _cluster_mmseqs2_fast."""
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="mmseqs-fast")
        cfg = pipeline.config

        mock_cluster = mocker.patch.object(
            pipeline, "_cluster_mmseqs2_fast",
            side_effect=lambda input_file, clusters_file: (
                clusters_file.write_text("#Sample\tNeighbors\tCluster\ns1\t2\t0\ns2\t2\t0\n"),
                clusters_file,
            )[-1],
        )

        def fake_makespan(**kw):
            out = Path(kw["output_file"])
            out.write_text("#Sample\tPartition\ns1\t1\ns2\t2\n")
            return Path(kw["stats_file"])
        mocker.patch.object(pipeline.makespan, "run_makespan", side_effect=fake_makespan)

        outputs = pipeline.run()
        mock_cluster.assert_called_once()
        assert "clusters" in outputs
        assert "makespan" in outputs


# ── _run_reduce pipeline ─────────────────────────────────────────────


class TestRunReduce:
    def test_reduce_pipeline(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, action="reduce")
        cfg = pipeline.config

        # Mock subprocess.run for KMA index (writes stdout to cluster file)
        mocker.patch("subprocess.run")

        def fake_makespan(**kw):
            out = Path(kw["output_file"])
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text("#Sample\tPartition\ns1\t1\n")
            stats = Path(kw["stats_file"])
            stats.write_text("stats\n")
            return stats
        mocker.patch.object(pipeline.makespan, "run_makespan", side_effect=fake_makespan)

        outputs = pipeline.run()
        assert "clusters" in outputs
        assert "makespan" in outputs

    def test_reduce_cleanup(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, action="reduce", keep_tmp=False)
        cfg = pipeline.config

        mocker.patch("subprocess.run")

        def fake_makespan(**kw):
            out = Path(kw["output_file"])
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text("x\n")
            return Path(kw["stats_file"])
        mocker.patch.object(pipeline.makespan, "run_makespan", side_effect=fake_makespan)

        mock_rmtree = mocker.patch("shutil.rmtree")
        pipeline.run()
        mock_rmtree.assert_called_once()


# ── Individual distance methods ──────────────────────────────────────


class TestDistanceKma:
    def test_kma_index_and_dist(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="cosine")
        mocker.patch.object(pipeline.kma_index, "run")
        mocker.patch.object(pipeline.kma_dist, "run")

        out = pipeline.config.tmp_dir / "out.phy"
        result = pipeline._distance_kma(
            pipeline.config.input_path, out,
        )
        assert result == out
        pipeline.kma_index.run.assert_called_once()
        pipeline.kma_dist.run.assert_called_once()


class TestDistanceMash:
    def test_mash_requires_max_length(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="mash", max_length=None)
        with pytest.raises(ValueError, match="max_length"):
            pipeline._distance_mash(Path("in.fsa"), Path("out.phy"))

    def test_mash_runs_triangle(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="mash", max_length=500)
        mocker.patch.object(pipeline.mash, "run_triangle")

        out = pipeline.config.tmp_dir / "out.phy"
        result = pipeline._distance_mash(pipeline.config.input_path, out)
        assert result == out
        pipeline.mash.run_triangle.assert_called_once()


class TestDistanceGgsearch:
    def test_ggsearch_requires_max_length(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="identity", max_length=None)
        with pytest.raises(ValueError, match="max_length"):
            pipeline._distance_ggsearch(Path("in.fsa"), Path("out.phy"))

    def test_ggsearch_builds_matrix(self, tmp_path, mocker):
        """Full ggsearch all-vs-all distance matrix construction."""
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="identity", max_length=500)
        cfg = pipeline.config

        # Mock run_alignment to create a result file
        def fake_alignment(**kw):
            out = Path(kw["output_file"])
            out.write_text(
                "The best scores are:\n"
                "s1   desc   (100)  150   0.90   100   50   20   100\n"
                ">>>\n"
            )
            return out
        mocker.patch.object(pipeline.ggsearch, "run_alignment", side_effect=fake_alignment)

        out = cfg.tmp_dir / "out.phy"
        result = pipeline._distance_ggsearch(cfg.input_path, out)
        assert result == out
        content = out.read_text()
        assert "\t2\n" in content  # 2 sequences
        assert "s1" in content
        assert "s2" in content


class TestDistanceMmseqs2:
    def test_mmseqs2_builds_matrix(self, tmp_path, mocker):
        """Full mmseqs2 all-vs-all distance matrix construction."""
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="mmseqs2")
        cfg = pipeline.config

        # Mock run_search to create a .m8 result
        def fake_search(**kw):
            m8 = Path(str(kw["out_prefix"]) + ".m8")
            m8.parent.mkdir(parents=True, exist_ok=True)
            m8.write_text(
                "s1\ts2\t0.90\t100\t10\t0\t1\t100\t1\t100\t1e-50\t200\n"
                "s2\ts1\t0.90\t100\t10\t0\t1\t100\t1\t100\t1e-50\t200\n"
            )
            return m8
        mocker.patch.object(pipeline.mmseqs2, "run_search", side_effect=fake_search)

        out = cfg.tmp_dir / "out.phy"
        result = pipeline._distance_mmseqs2(cfg.input_path, out)
        assert result == out
        content = out.read_text()
        assert "\t2\n" in content  # 2 sequences
        assert "s1" in content
        assert "s2" in content
        # Distance = 1 - 0.90 ≈ 0.1
        assert "s2\t" in content  # s2 has a distance row


# ── MMseqs-fast clustering ────────────────────────────────────────────


class TestClusterMmseqsFast:
    def test_cluster_and_convert(self, tmp_path, mocker):
        """_cluster_mmseqs2_fast runs easy-cluster and converts to DBSCAN format."""
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="mmseqs-fast")
        cfg = pipeline.config

        # Mock run_cluster to create a cluster TSV
        def fake_run_cluster(**kw):
            tsv = Path(str(kw["out_prefix"]) + "_cluster.tsv")
            tsv.parent.mkdir(parents=True, exist_ok=True)
            tsv.write_text("rep1\ts1\nrep1\ts2\nrep1\trep1\nrep2\ts3\nrep2\trep2\n")
            return tsv
        mocker.patch.object(
            pipeline.mmseqs2_cluster, "run_cluster", side_effect=fake_run_cluster
        )

        clusters_file = cfg.results_dir / "clusters.tsv"
        result = pipeline._cluster_mmseqs2_fast(cfg.input_path, clusters_file)
        assert result == clusters_file

        content = clusters_file.read_text()
        assert "#Sample\tNeighbors\tCluster" in content
        lines = content.strip().split("\n")
        assert len(lines) == 6  # header + 5 members
        # Cluster 0 has 3 members, cluster 1 has 2
        assert "s1\t3\t0" in content
        assert "s2\t3\t0" in content
        assert "rep1\t3\t0" in content
        assert "s3\t2\t1" in content
        assert "rep2\t2\t1" in content


# ── Hobohm1 reduction ────────────────────────────────────────────────


class TestRunHobohm1:
    def test_hobohm1_cdhit(self, tmp_path, mocker):
        pipeline = _make_pipeline(
            tmp_path, mocker,
            approach="hobohm_reduce",
            hobohm1_distance=0.05,
            hobohm1_method="cdhit",
        )
        mocker.patch.object(pipeline.cdhit, "run")

        result = pipeline._run_hobohm1(
            input_file=pipeline.config.input_path,
            output_prefix=pipeline.config.tmp_dir / "hobohm",
        )
        assert str(result).endswith(".fsa")
        pipeline.cdhit.run.assert_called_once()

    def test_hobohm1_kma(self, tmp_path, mocker):
        pipeline = _make_pipeline(
            tmp_path, mocker,
            approach="hobohm_split",
            hobohm1_distance=0.05,
            hobohm1_method="kma",
        )
        mocker.patch("subprocess.run")

        result = pipeline._run_hobohm1(
            input_file=pipeline.config.input_path,
            output_prefix=pipeline.config.tmp_dir / "hobohm",
        )
        # KMA hobohm returns the original input
        assert result == pipeline.config.input_path


# ── Merge tables ─────────────────────────────────────────────────────


class TestMergeTables:
    def test_merge(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker)

        clusters = tmp_path / "clusters.tsv"
        clusters.write_text("#Sample\tCluster\nseq1\t1\nseq2\t2\n")

        makespan = tmp_path / "makespan.tsv"
        makespan.write_text("#Sample\tPartition\nseq1\t1\nseq2\t2\n")

        out = tmp_path / "merged.tsv"
        result = pipeline._merge_tables(clusters, makespan, out)
        assert result.exists()
        content = result.read_text()
        assert "seq1" in content
        assert "Cluster" in content
        assert "Partition" in content


# ── Create FASTA partitions ──────────────────────────────────────────


class TestCreateFastaPartitions:
    def test_split_fasta(self, tmp_path, mocker):
        fasta = tmp_path / "input.fsa"
        cfg = _make_config(tmp_path, input_path=fasta, bins=2)
        # Write 3 sequences AFTER _make_config (which overwrites input.fsa)
        fasta.write_text(">s1\nATCG\n>s2\nGGCC\n>s3\nTTTT\n")

        fake_exe = tmp_path / "tool"
        fake_exe.write_text("#!/bin/sh\n")
        fake_exe.chmod(0o755)
        mocker.patch.object(
            SpanSeqPipeline, "_resolve_tool", return_value=fake_exe
        )
        pipeline = SpanSeqPipeline(cfg)

        partitions = tmp_path / "partitions.tsv"
        partitions.write_text("id\tpartition\ns1\t1\ns2\t2\ns3\t1\n")

        result = pipeline._create_fasta_partitions(partitions, fasta)
        assert len(result) == 2
        p1_content = result[0].read_text()
        assert "s1" in p1_content
        assert "s3" in p1_content
        p2_content = result[1].read_text()
        assert "s2" in p2_content

    def test_split_fasta_description_fallback(self, tmp_path, mocker):
        """When id doesn't match, falls back to full description."""
        fasta = tmp_path / "input.fsa"
        cfg = _make_config(tmp_path, input_path=fasta, bins=2)
        fasta.write_text(">gene_A some desc\nATCG\n>gene_B other\nGGCC\n")

        fake_exe = tmp_path / "tool"
        fake_exe.write_text("#!/bin/sh\n")
        fake_exe.chmod(0o755)
        mocker.patch.object(
            SpanSeqPipeline, "_resolve_tool", return_value=fake_exe
        )
        pipeline = SpanSeqPipeline(cfg)

        partitions = tmp_path / "partitions.tsv"
        # Use full description as id
        partitions.write_text(
            "id\tpartition\ngene_A some desc\t1\ngene_B other\t2\n"
        )

        result = pipeline._create_fasta_partitions(partitions, fasta)
        p1 = result[0].read_text()
        assert "gene_A" in p1
        p2 = result[1].read_text()
        assert "gene_B" in p2

    def test_split_fasta_missing_sequence_warning(self, tmp_path, mocker):
        """Sequences not in partition table are skipped with warning."""
        fasta = tmp_path / "input.fsa"
        cfg = _make_config(tmp_path, input_path=fasta, bins=2)
        fasta.write_text(">s1\nATCG\n>unknown\nGGCC\n")

        fake_exe = tmp_path / "tool"
        fake_exe.write_text("#!/bin/sh\n")
        fake_exe.chmod(0o755)
        mocker.patch.object(
            SpanSeqPipeline, "_resolve_tool", return_value=fake_exe
        )
        pipeline = SpanSeqPipeline(cfg)

        partitions = tmp_path / "partitions.tsv"
        partitions.write_text("id\tpartition\ns1\t1\n")

        result = pipeline._create_fasta_partitions(partitions, fasta)
        p1 = result[0].read_text()
        assert "s1" in p1
        # unknown is skipped
        p2 = result[1].read_text()
        assert p2 == ""


# ── Add class columns ────────────────────────────────────────────────


class TestAddClassColumns:
    def test_add_class_columns(self, tmp_path, mocker):
        imb = tmp_path / "imb.tsv"
        imb.write_text("s1\tA\ns2\tB\n")
        pipeline = _make_pipeline(tmp_path, mocker, imbalance_file=imb)

        clusters = tmp_path / "clusters.tsv"
        clusters.write_text("#Sample\tCluster\ns1\t1\ns2\t2\n")

        out = tmp_path / "clusters_class.tsv"
        result = pipeline._add_class_columns(clusters, out)
        assert result.exists()
        content = result.read_text()
        assert "Class_1" in content

    def test_add_class_columns_fewer_raises(self, tmp_path, mocker):
        """Imbalance file with fewer sequences raises ValueError."""
        imb = tmp_path / "imb.tsv"
        imb.write_text("s1\tA\n")  # only s1, no s2
        pipeline = _make_pipeline(tmp_path, mocker, imbalance_file=imb)

        clusters = tmp_path / "clusters.tsv"
        clusters.write_text("#Sample\tCluster\ns1\t1\ns2\t2\n")

        out = tmp_path / "clusters_class.tsv"
        with pytest.raises(ValueError, match="fewer sequences"):
            pipeline._add_class_columns(clusters, out)


# ── Validation ───────────────────────────────────────────────────────


class TestDistanceValidation:
    def test_mash_requires_max_length(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="mash", max_length=None)
        with pytest.raises(ValueError, match="max_length"):
            pipeline._distance_mash(Path("in.fsa"), Path("out.phy"))

    def test_ggsearch_requires_max_length(self, tmp_path, mocker):
        pipeline = _make_pipeline(tmp_path, mocker, distance_method="identity", max_length=None)
        with pytest.raises(ValueError, match="max_length"):
            pipeline._distance_ggsearch(Path("in.fsa"), Path("out.phy"))
