"""Tests for GGSearchApp (AlignerRunner)."""
import pytest
from pathlib import Path
from cgecore.applications.base import AlignerRunner
from spanseq.applications.ggsearch import GGSearchApp


class TestGGSearchApp:
    def test_is_aligner_runner(self):
        assert issubclass(GGSearchApp, AlignerRunner)

    def test_nucleotide_command(self, fake_exec, tmp_path):
        app = GGSearchApp(exec_path=fake_exec)

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
            seq_type="nucleotides", max_len=500, cores=2,
        )
        assert str(fake_exec) in cmd[0]
        assert "-n" in cmd
        assert "-T" in cmd
        idx = cmd.index("-T")
        assert cmd[idx + 1] == "2"

    def test_aminoacid_flag(self, fake_exec, tmp_path):
        app = GGSearchApp(exec_path=fake_exec)

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
        assert "-p" in cmd
        assert "-n" not in cmd

    def test_invalid_seq_type(self, fake_exec, tmp_path):
        app = GGSearchApp(exec_path=fake_exec)

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

    def test_map_outputs(self, fake_exec):
        app = GGSearchApp(exec_path=fake_exec)
        outputs = app.map_outputs(workdir=Path("/tmp"), out_prefix="/tmp/out")
        assert outputs["alignment"] == Path("/tmp/out.ggsearch36")

    def test_map_outputs_none(self, fake_exec):
        app = GGSearchApp(exec_path=fake_exec)
        assert app.map_outputs(workdir=Path("/tmp"), out_prefix=None) == {}

    def test_run_alignment(self, fake_exec, tmp_path, mocker):
        app = GGSearchApp(exec_path=fake_exec)
        query = tmp_path / "query.fsa"
        query.write_text(">q\nATCG\n")
        target = tmp_path / "target.fsa"
        target.write_text(">t\nATCG\n")

        mocker.patch("subprocess.run")
        out = tmp_path / "out.ggsearch36"
        result = app.run_alignment(
            output_file=out,
            query_file=query,
            target_file=target,
            seq_type="nucleotides",
            max_len=500,
        )
        assert result == out


class TestParseResultFile:
    def test_parse_sample_output(self, tmp_path):
        output = tmp_path / "result.ggsearch36"
        output.write_text(
            "some header lines\n"
            "The best scores are:\n"
            "gene_A   description   (100)  150   0.95   100   50   20   100\n"
            "gene_B   description   (200)  120   0.80   80    40   10   200\n"
            ">>>\n"
        )
        results = GGSearchApp.parse_result_file(output)
        assert "gene_A" in results
        assert "gene_B" in results

    def test_parse_missing_file(self, tmp_path):
        result = GGSearchApp.parse_result_file(tmp_path / "nonexistent")
        assert result == {}

    def test_parse_empty_file(self, tmp_path):
        empty = tmp_path / "empty.ggsearch36"
        empty.write_text("")
        result = GGSearchApp.parse_result_file(empty)
        assert result == {}

    def test_parse_short_score_line(self, tmp_path):
        """Lines with fewer than 6 fields are skipped."""
        output = tmp_path / "result.ggsearch36"
        output.write_text(
            "The best scores are:\n"
            "gene_A   50   0.95\n"
            ">>>\n"
        )
        result = GGSearchApp.parse_result_file(output)
        assert result == {}

    def test_parse_empty_line_in_scores(self, tmp_path):
        """Blank lines in the scores section are skipped."""
        output = tmp_path / "result.ggsearch36"
        output.write_text(
            "The best scores are:\n"
            "\n"
            "gene_A   description   (100)  150   0.95   100   50   20   100\n"
            ">>>\n"
        )
        result = GGSearchApp.parse_result_file(output)
        assert "gene_A" in result

    def test_parse_malformed_values(self, tmp_path):
        """Non-numeric score values are skipped."""
        output = tmp_path / "result.ggsearch36"
        output.write_text(
            "The best scores are:\n"
            "gene_A   desc   (100)  XXX   YYY   ZZZ   AAA   BBB   CCC\n"
            ">>>\n"
        )
        result = GGSearchApp.parse_result_file(output)
        assert result == {}
