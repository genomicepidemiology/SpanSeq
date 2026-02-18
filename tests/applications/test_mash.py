"""Tests for MashApp command building."""
from pathlib import Path
from spanseq.applications.mash import MashApp


class TestMashApp:
    def test_nucleotide_command(self, fake_exec):
        app = MashApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/data/seqs.fsa",
            seq_type="nucleotides",
            kmer_size=7,
            sketch_size=1000,
            threads=4,
        )
        assert cmd[0] == str(fake_exec)
        assert "triangle" in cmd
        assert "-a" not in cmd
        idx = cmd.index("-k")
        assert cmd[idx + 1] == "7"
        idx = cmd.index("-s")
        assert cmd[idx + 1] == "1000"
        idx = cmd.index("-p")
        assert cmd[idx + 1] == "4"

    def test_aminoacid_flag(self, fake_exec):
        app = MashApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/data/seqs.fsa",
            seq_type="aminoacids",
        )
        assert "-a" in cmd

    def test_batch_input(self, fake_exec):
        app = MashApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/data/batch.txt",
            input_format="batch",
        )
        assert "-l" in cmd

    def test_extra_args_sketch(self, fake_exec):
        app = MashApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/data/seqs.fsa",
            extra_args_sketch="-S 42",
        )
        assert "-S" in cmd
        assert "42" in cmd

    def test_extra_args_triangle(self, fake_exec):
        app = MashApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/data/seqs.fsa",
            extra_args_triangle="-E 0.01",
        )
        assert "-E" in cmd
        assert "0.01" in cmd

    def test_run_triangle(self, fake_exec, tmp_path, mocker):
        app = MashApp(exec_path=fake_exec)
        mocker.patch("subprocess.run")
        out = tmp_path / "out.phy"
        result = app.run_triangle(
            output_file=out,
            input_file="/data/seqs.fsa",
        )
        assert result == out

    def test_map_outputs(self, fake_exec):
        app = MashApp(exec_path=fake_exec)
        outputs = app.map_outputs(
            workdir=Path("/tmp"),
            out_prefix="/tmp/out",
        )
        assert outputs["phy"] == Path("/tmp/out.phy")

    def test_map_outputs_none(self, fake_exec):
        app = MashApp(exec_path=fake_exec)
        assert app.map_outputs(workdir=Path("/tmp"), out_prefix=None) == {}
