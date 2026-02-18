"""Tests for KmaIndexApp and KmaDistApp command building."""
from pathlib import Path
from spanseq.applications.kma import KmaIndexApp, KmaDistApp


class TestKmaIndexApp:
    def test_basic_command(self, fake_exec):
        app = KmaIndexApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/data/seqs.fsa",
            output_prefix="/tmp/index",
        )
        assert cmd[0] == str(fake_exec)
        assert "index" in cmd
        assert "-i" in cmd
        assert "/data/seqs.fsa" in cmd
        assert "-o" in cmd
        assert "/tmp/index" in cmd

    def test_batch_input(self, fake_exec):
        app = KmaIndexApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/data/batch.txt",
            output_prefix="/tmp/index",
            input_format="batch",
        )
        assert "-batch" in cmd
        assert "-i" not in cmd

    def test_kmer_size(self, fake_exec):
        app = KmaIndexApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/data/seqs.fsa",
            output_prefix="/tmp/index",
            kmer_size=16,
        )
        idx = cmd.index("-k")
        assert cmd[idx + 1] == "16"

    def test_hobohm1_flags(self, fake_exec):
        app = KmaIndexApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/data/seqs.fsa",
            output_prefix="/tmp/index",
            hobohm1_value=0.8,
        )
        assert "-hq" in cmd
        assert "-ht" in cmd
        assert "-and" in cmd
        idx = cmd.index("-hq")
        assert cmd[idx + 1] == "0.8"

    def test_megadb_flag(self, fake_exec):
        app = KmaIndexApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/data/seqs.fsa",
            output_prefix="/tmp/index",
            megadb=True,
        )
        assert "-ME" in cmd

    def test_map_outputs(self, fake_exec):
        app = KmaIndexApp(exec_path=fake_exec)
        outputs = app.map_outputs(
            workdir=Path("/tmp"),
            out_prefix="/tmp/index",
        )
        assert "comp" in outputs
        assert "length" in outputs
        assert "name" in outputs
        assert "seq" in outputs


class TestKmaDistApp:
    def test_basic_command(self, fake_exec):
        app = KmaDistApp(exec_path=fake_exec)
        cmd = app.build_command(
            db_prefix="/tmp/index",
            output_file="/tmp/out.phy",
            method=256,
            threads=4,
        )
        assert cmd[0] == str(fake_exec)
        assert "dist" in cmd
        assert "-d" in cmd
        idx = cmd.index("-d")
        assert cmd[idx + 1] == "256"
        idx = cmd.index("-t")
        assert cmd[idx + 1] == "4"

    def test_tmp_dir(self, fake_exec):
        app = KmaDistApp(exec_path=fake_exec)
        cmd = app.build_command(
            db_prefix="/tmp/index",
            output_file="/tmp/out.phy",
            tmp_dir="/tmp/scratch",
        )
        assert "-tmp" in cmd
        assert "/tmp/scratch" in cmd

    def test_map_outputs(self, fake_exec):
        app = KmaDistApp(exec_path=fake_exec)
        outputs = app.map_outputs(
            workdir=Path("/tmp"),
            out_prefix="/tmp/out",
        )
        assert outputs["phy"] == Path("/tmp/out.phy")
