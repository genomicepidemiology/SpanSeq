"""Tests for CCPhyloDbscanApp, CCPhyloMakespanApp, and CCPhyloTreeApp."""
from pathlib import Path
from spanseq.applications.ccphylo import CCPhyloDbscanApp, CCPhyloMakespanApp, CCPhyloTreeApp


class TestDbscanApp:
    def test_basic_command(self, fake_exec):
        app = CCPhyloDbscanApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/tmp/dist.phy",
            output_file="/tmp/clusters.tsv",
            dist_value=0.3,
        )
        assert cmd[0] == str(fake_exec)
        assert "dbscan" in cmd
        assert "-p" in cmd
        idx = cmd.index("-e")
        assert cmd[idx + 1] == "0.3"

    def test_memory_disk(self, fake_exec):
        app = CCPhyloDbscanApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/tmp/dist.phy",
            output_file="/tmp/clusters.tsv",
            dist_value=0.3,
            memory_disk=True,
        )
        assert "-H" in cmd

    def test_tmp_dir(self, fake_exec):
        app = CCPhyloDbscanApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/tmp/dist.phy",
            output_file="/tmp/clusters.tsv",
            dist_value=0.3,
            tmp_dir="/scratch",
        )
        assert "-T" in cmd
        assert "/scratch" in cmd

    def test_extra_args(self, fake_exec):
        app = CCPhyloDbscanApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/tmp/dist.phy",
            output_file="/tmp/clusters.tsv",
            dist_value=0.3,
            extra_args="-v --debug",
        )
        assert "-v" in cmd
        assert "--debug" in cmd

    def test_map_outputs(self, fake_exec):
        app = CCPhyloDbscanApp(exec_path=fake_exec)
        outputs = app.map_outputs(workdir=Path("/tmp"), out_prefix="/tmp/out")
        assert outputs["clusters"] == Path("/tmp/out_clusters.tsv")

    def test_map_outputs_none(self, fake_exec):
        app = CCPhyloDbscanApp(exec_path=fake_exec)
        assert app.map_outputs(workdir=Path("/tmp"), out_prefix=None) == {}


class TestMakespanApp:
    def test_basic_command(self, fake_exec):
        app = CCPhyloMakespanApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/tmp/clusters.tsv",
            output_file="/tmp/makespan.tsv",
            stats_file="/tmp/stats.tsv",
            machines=5,
        )
        assert cmd[0] == str(fake_exec)
        assert "makespan" in cmd
        idx = cmd.index("-l")
        assert cmd[idx + 1] == "5"

    def test_method_and_weights(self, fake_exec):
        app = CCPhyloMakespanApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/tmp/clusters.tsv",
            output_file="/tmp/makespan.tsv",
            stats_file="/tmp/stats.tsv",
            machines=3,
            method="DFF",
            weight_method="log2",
        )
        idx = cmd.index("-m")
        assert cmd[idx + 1] == "DFF"
        idx = cmd.index("-w")
        assert cmd[idx + 1] == "log2"

    def test_class_columns(self, fake_exec):
        app = CCPhyloMakespanApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/tmp/clusters.tsv",
            output_file="/tmp/makespan.tsv",
            stats_file="/tmp/stats.tsv",
            machines=3,
            class_columns="4,5",
        )
        assert "-c" in cmd
        assert "4,5" in cmd

    def test_extra_args(self, fake_exec):
        app = CCPhyloMakespanApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/tmp/clusters.tsv",
            output_file="/tmp/makespan.tsv",
            stats_file="/tmp/stats.tsv",
            machines=3,
            extra_args="--verbose",
        )
        assert "--verbose" in cmd

    def test_run_makespan(self, fake_exec, tmp_path, mocker):
        app = CCPhyloMakespanApp(exec_path=fake_exec)
        mocker.patch("subprocess.run")
        stats = tmp_path / "stats.tsv"
        result = app.run_makespan(
            stats_file=stats,
            input_file="/tmp/clusters.tsv",
            output_file="/tmp/makespan.tsv",
            machines=3,
        )
        assert result == stats

    def test_map_outputs(self, fake_exec):
        app = CCPhyloMakespanApp(exec_path=fake_exec)
        outputs = app.map_outputs(workdir=Path("/tmp"), out_prefix="/tmp/out")
        assert outputs["makespan"] == Path("/tmp/out_makespan.tsv")
        assert outputs["stats"] == Path("/tmp/out_makespan_stats.tsv")

    def test_map_outputs_none(self, fake_exec):
        app = CCPhyloMakespanApp(exec_path=fake_exec)
        assert app.map_outputs(workdir=Path("/tmp"), out_prefix=None) == {}


class TestTreeApp:
    def test_basic_command(self, fake_exec):
        app = CCPhyloTreeApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/tmp/dist.phy",
            output_file="/tmp/tree.nwk",
        )
        assert cmd[0] == str(fake_exec)
        assert "tree" in cmd
        idx = cmd.index("-m")
        assert cmd[idx + 1] == "dnj"
        idx = cmd.index("-t")
        assert cmd[idx + 1] == "1"

    def test_method_and_threads(self, fake_exec):
        app = CCPhyloTreeApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/tmp/dist.phy",
            output_file="/tmp/tree.nwk",
            method="nj",
            threads=8,
        )
        idx = cmd.index("-m")
        assert cmd[idx + 1] == "nj"
        idx = cmd.index("-t")
        assert cmd[idx + 1] == "8"

    def test_memory_disk(self, fake_exec):
        app = CCPhyloTreeApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/tmp/dist.phy",
            output_file="/tmp/tree.nwk",
            memory_disk=True,
        )
        assert "-H" in cmd

    def test_tmp_dir(self, fake_exec):
        app = CCPhyloTreeApp(exec_path=fake_exec)
        cmd = app.build_command(
            input_file="/tmp/dist.phy",
            output_file="/tmp/tree.nwk",
            tmp_dir="/scratch",
        )
        assert "-T" in cmd
        assert "/scratch" in cmd

    def test_map_outputs(self, fake_exec):
        app = CCPhyloTreeApp(exec_path=fake_exec)
        outputs = app.map_outputs(workdir=Path("/tmp"), out_prefix="/tmp/out")
        assert outputs["tree"] == Path("/tmp/out.nwk")

    def test_map_outputs_none(self, fake_exec):
        app = CCPhyloTreeApp(exec_path=fake_exec)
        assert app.map_outputs(workdir=Path("/tmp"), out_prefix=None) == {}
