"""Tests for MMseqs2SearchApp and MMseqs2ClusterApp."""
import pytest
from pathlib import Path
from cgecore.applications.base import ApplicationRunner, AlignerRunner
from spanseq.applications.mmseqs2 import MMseqs2SearchApp, MMseqs2ClusterApp


class TestMMseqs2SearchApp:
    def test_is_aligner_runner(self):
        assert issubclass(MMseqs2SearchApp, AlignerRunner)

    def test_nucleotide_command(self, fake_exec, tmp_path):
        app = MMseqs2SearchApp(exec_path=fake_exec)

        query = tmp_path / "query.fsa"
        query.write_text(">q\nATCG\n")
        target = tmp_path / "target.fsa"
        target.write_text(">t\nATCG\n")

        from cgecore.signatures.sample import SampleSpec
        from cgecore.signatures.database import DatabaseSpec

        q = SampleSpec(type="assembled", files=[query])
        db = DatabaseSpec(
            name=target.name, root=target.parent,
            indexed=False, indexer=None, sequencetype="genes",
        )

        cmd = app.build_command(
            query=q, db=db, out_prefix=str(tmp_path / "out"),
            seq_type="nucleotides", cores=4,
        )
        assert "easy-search" in cmd
        assert "--search-type" in cmd
        idx = cmd.index("--search-type")
        assert cmd[idx + 1] == "3"
        idx = cmd.index("--threads")
        assert cmd[idx + 1] == "4"

    def test_aminoacid_search_type(self, fake_exec, tmp_path):
        app = MMseqs2SearchApp(exec_path=fake_exec)

        query = tmp_path / "query.fsa"
        query.write_text(">q\nMLLK\n")
        target = tmp_path / "target.fsa"
        target.write_text(">t\nMLLK\n")

        from cgecore.signatures.sample import SampleSpec
        from cgecore.signatures.database import DatabaseSpec

        q = SampleSpec(type="assembled", files=[query])
        db = DatabaseSpec(
            name=target.name, root=target.parent,
            indexed=False, indexer=None, sequencetype="genes",
        )

        cmd = app.build_command(
            query=q, db=db, out_prefix=str(tmp_path / "out"),
            seq_type="aminoacids",
        )
        idx = cmd.index("--search-type")
        assert cmd[idx + 1] == "1"

    def test_invalid_seq_type(self, fake_exec, tmp_path):
        app = MMseqs2SearchApp(exec_path=fake_exec)

        query = tmp_path / "query.fsa"
        query.write_text(">q\nA\n")
        target = tmp_path / "target.fsa"
        target.write_text(">t\nA\n")

        from cgecore.signatures.sample import SampleSpec
        from cgecore.signatures.database import DatabaseSpec

        q = SampleSpec(type="assembled", files=[query])
        db = DatabaseSpec(
            name=target.name, root=target.parent,
            indexed=False, indexer=None, sequencetype="genes",
        )

        with pytest.raises(ValueError, match="seq_type"):
            app.build_command(
                query=q, db=db, out_prefix=str(tmp_path / "out"),
                seq_type="rna",
            )

    def test_min_seq_id_omitted_when_zero(self, fake_exec, tmp_path):
        app = MMseqs2SearchApp(exec_path=fake_exec)

        query = tmp_path / "query.fsa"
        query.write_text(">q\nA\n")
        target = tmp_path / "target.fsa"
        target.write_text(">t\nA\n")

        from cgecore.signatures.sample import SampleSpec
        from cgecore.signatures.database import DatabaseSpec

        q = SampleSpec(type="assembled", files=[query])
        db = DatabaseSpec(
            name=target.name, root=target.parent,
            indexed=False, indexer=None, sequencetype="genes",
        )

        cmd = app.build_command(
            query=q, db=db, out_prefix=str(tmp_path / "out"),
            min_seq_id=0.0,
        )
        assert "--min-seq-id" not in cmd

    def test_min_seq_id_added_when_positive(self, fake_exec, tmp_path):
        app = MMseqs2SearchApp(exec_path=fake_exec)

        query = tmp_path / "query.fsa"
        query.write_text(">q\nA\n")
        target = tmp_path / "target.fsa"
        target.write_text(">t\nA\n")

        from cgecore.signatures.sample import SampleSpec
        from cgecore.signatures.database import DatabaseSpec

        q = SampleSpec(type="assembled", files=[query])
        db = DatabaseSpec(
            name=target.name, root=target.parent,
            indexed=False, indexer=None, sequencetype="genes",
        )

        cmd = app.build_command(
            query=q, db=db, out_prefix=str(tmp_path / "out"),
            min_seq_id=0.3,
        )
        assert "--min-seq-id" in cmd
        idx = cmd.index("--min-seq-id")
        assert cmd[idx + 1] == "0.3"

    def test_map_outputs(self, fake_exec):
        app = MMseqs2SearchApp(exec_path=fake_exec)
        outputs = app.map_outputs(workdir=Path("/tmp"), out_prefix="/tmp/out")
        assert outputs["result"] == Path("/tmp/out.m8")

    def test_map_outputs_none(self, fake_exec):
        app = MMseqs2SearchApp(exec_path=fake_exec)
        assert app.map_outputs(workdir=Path("/tmp"), out_prefix=None) == {}


class TestParseResultFile:
    def test_parse_m8_output(self, tmp_path):
        m8 = tmp_path / "result.m8"
        m8.write_text(
            "seq1\tseq2\t0.95\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200\n"
            "seq1\tseq3\t0.80\t90\t18\t0\t1\t90\t1\t90\t1e-30\t150\n"
            "seq2\tseq3\t0.70\t80\t24\t0\t1\t80\t1\t80\t1e-20\t100\n"
        )
        results = MMseqs2SearchApp.parse_result_file(m8)
        assert results["seq1"]["seq2"]["identity"] == 0.95
        assert results["seq1"]["seq3"]["identity"] == 0.80
        assert results["seq2"]["seq3"]["identity"] == 0.70

    def test_keeps_highest_identity(self, tmp_path):
        m8 = tmp_path / "result.m8"
        m8.write_text(
            "seq1\tseq2\t0.80\t50\t10\t0\t1\t50\t1\t50\t1e-10\t80\n"
            "seq1\tseq2\t0.95\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200\n"
        )
        results = MMseqs2SearchApp.parse_result_file(m8)
        assert results["seq1"]["seq2"]["identity"] == 0.95

    def test_parse_missing_file(self, tmp_path):
        result = MMseqs2SearchApp.parse_result_file(tmp_path / "nope.m8")
        assert result == {}

    def test_parse_empty_file(self, tmp_path):
        m8 = tmp_path / "empty.m8"
        m8.write_text("")
        result = MMseqs2SearchApp.parse_result_file(m8)
        assert result == {}

    def test_parse_comment_lines_skipped(self, tmp_path):
        m8 = tmp_path / "result.m8"
        m8.write_text(
            "# header comment\n"
            "seq1\tseq2\t0.95\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200\n"
        )
        results = MMseqs2SearchApp.parse_result_file(m8)
        assert results["seq1"]["seq2"]["identity"] == 0.95

    def test_parse_short_line_skipped(self, tmp_path):
        m8 = tmp_path / "result.m8"
        m8.write_text("seq1\tseq2\n")
        result = MMseqs2SearchApp.parse_result_file(m8)
        assert result == {}

    def test_parse_malformed_identity(self, tmp_path):
        m8 = tmp_path / "result.m8"
        m8.write_text("seq1\tseq2\tNOT_A_NUMBER\t100\n")
        result = MMseqs2SearchApp.parse_result_file(m8)
        assert result == {}


class TestMMseqs2ClusterApp:
    def test_is_application_runner(self):
        assert issubclass(MMseqs2ClusterApp, ApplicationRunner)

    def test_nucleotide_command(self, fake_exec, tmp_path):
        app = MMseqs2ClusterApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file=tmp_path / "input.fsa",
            out_prefix=tmp_path / "out",
            seq_type="nucleotides",
            min_seq_id=0.7,
            threads=4,
        )
        assert "easy-cluster" in cmd
        assert "--search-type" in cmd
        idx = cmd.index("--search-type")
        assert cmd[idx + 1] == "3"
        idx = cmd.index("--min-seq-id")
        assert cmd[idx + 1] == "0.7"
        idx = cmd.index("--threads")
        assert cmd[idx + 1] == "4"

    def test_aminoacid_command(self, fake_exec, tmp_path):
        app = MMseqs2ClusterApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file=tmp_path / "input.fsa",
            out_prefix=tmp_path / "out",
            seq_type="aminoacids",
        )
        idx = cmd.index("--search-type")
        assert cmd[idx + 1] == "1"

    def test_invalid_seq_type(self, fake_exec, tmp_path):
        app = MMseqs2ClusterApp(exec_path=fake_exec)
        with pytest.raises(ValueError, match="seq_type"):
            app.build_command(
                input_file=tmp_path / "input.fsa",
                out_prefix=tmp_path / "out",
                seq_type="rna",
            )

    def test_map_outputs(self, fake_exec):
        app = MMseqs2ClusterApp(exec_path=fake_exec)
        outputs = app.map_outputs(workdir=Path("/tmp"), out_prefix="/tmp/out")
        assert outputs["cluster_tsv"] == Path("/tmp/out_cluster.tsv")

    def test_map_outputs_none(self, fake_exec):
        app = MMseqs2ClusterApp(exec_path=fake_exec)
        assert app.map_outputs(workdir=Path("/tmp"), out_prefix=None) == {}

    def test_run_cluster(self, fake_exec, tmp_path, mocker):
        app = MMseqs2ClusterApp(exec_path=fake_exec)
        mocker.patch.object(app, "run")

        result = app.run_cluster(
            input_file=tmp_path / "input.fsa",
            out_prefix=tmp_path / "out",
            seq_type="nucleotides",
            min_seq_id=0.7,
            threads=2,
        )
        assert result == tmp_path / "out_cluster.tsv"
        app.run.assert_called_once()


class TestParseClusterTsv:
    def test_parse_basic(self, tmp_path):
        tsv = tmp_path / "cluster.tsv"
        tsv.write_text(
            "rep1\tmember1\n"
            "rep1\tmember2\n"
            "rep1\trep1\n"
            "rep2\tmember3\n"
            "rep2\trep2\n"
        )
        result = MMseqs2ClusterApp.parse_cluster_tsv(tsv)
        assert len(result) == 2
        assert set(result["rep1"]) == {"member1", "member2", "rep1"}
        assert set(result["rep2"]) == {"member3", "rep2"}

    def test_parse_missing_file(self, tmp_path):
        result = MMseqs2ClusterApp.parse_cluster_tsv(tmp_path / "nope.tsv")
        assert result == {}

    def test_parse_empty_file(self, tmp_path):
        tsv = tmp_path / "empty.tsv"
        tsv.write_text("")
        result = MMseqs2ClusterApp.parse_cluster_tsv(tsv)
        assert result == {}

    def test_parse_short_line_skipped(self, tmp_path):
        tsv = tmp_path / "cluster.tsv"
        tsv.write_text("rep1\n")
        result = MMseqs2ClusterApp.parse_cluster_tsv(tsv)
        assert result == {}

    def test_parse_single_cluster(self, tmp_path):
        tsv = tmp_path / "cluster.tsv"
        tsv.write_text("rep1\trep1\n")
        result = MMseqs2ClusterApp.parse_cluster_tsv(tsv)
        assert result == {"rep1": ["rep1"]}
