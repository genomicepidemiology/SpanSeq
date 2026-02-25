"""CCPhylo application runners for clustering (dbscan), partitioning (makespan), and tree building."""
from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional, List, Any
from cgecore.applications.base import ApplicationRunner

import logging

logger = logging.getLogger(__name__)


class CCPhyloDbscanApp(ApplicationRunner):
    """Runner for `ccphylo dbscan` — density-based clustering from a distance matrix."""

    def __init__(self, exec_path: Path | str = "ccphylo") -> None:
        super().__init__(exec_path=Path(exec_path), tool_name="ccphylo")

    def build_command(
        self,
        *,
        input_file: Path | str,
        output_file: Path | str,
        dist_value: float,
        memory_disk: bool = False,
        tmp_dir: Optional[Path | str] = None,
        extra_args: str = "",
        **kwargs: Any,
    ) -> List[str]:
        """Build `ccphylo dbscan` command.

        Args:
            input_file: Distance matrix (.phy).
            output_file: Output cluster file (.tsv).
            dist_value: Maximum distance threshold (-e).
            memory_disk: Allocate distance matrix on disk (-H).
            tmp_dir: Temporary directory (-T).
            extra_args: Additional CLI arguments.
        """
        cmd = [
            str(self.exec_path), "dbscan",
            "-i", str(input_file),
            "-o", str(output_file),
            "-e", str(dist_value),
            "-p",
        ]

        if memory_disk:
            cmd.append("-H")

        if tmp_dir is not None:
            cmd += ["-T", str(tmp_dir)]

        if extra_args:
            cmd += extra_args.split()

        return cmd

    def map_outputs(
        self,
        workdir: Path,
        *,
        out_prefix: Optional[str] = None,
        app_args: Optional[Dict[str, Any]] = None,
        **kwargs: Any,
    ) -> Dict[str, Path]:
        if out_prefix is None:
            return {}
        return {"clusters": Path(str(out_prefix) + "_clusters.tsv")}


class CCPhyloMakespanApp(ApplicationRunner):
    """Runner for `ccphylo makespan` — partition clusters into bins."""

    def __init__(self, exec_path: Path | str = "ccphylo") -> None:
        super().__init__(exec_path=Path(exec_path), tool_name="ccphylo")

    def build_command(
        self,
        *,
        input_file: Path | str,
        output_file: Path | str,
        stats_file: Path | str,
        machines: int | str,
        field_cluster: int = 3,
        method: str = "DBF",
        weight_method: str = "none",
        class_columns: Optional[str] = None,
        extra_args: str = "",
        **kwargs: Any,
    ) -> List[str]:
        """Build `ccphylo makespan` command.

        Args:
            input_file: Cluster file (.tsv from dbscan).
            output_file: Output partition file (.tsv).
            stats_file: Output statistics file (.tsv).
            machines: Number of bins (-l). Int or comma-separated string.
            field_cluster: Column index for cluster field (-k).
            method: Makespan method (-m): "DBF" or "DFF".
            weight_method: Weighting method (-w): "none", "logX", "powX", "expX".
            class_columns: Class columns for imbalance-aware partitioning (-c).
            extra_args: Additional CLI arguments.
        """
        # ccphylo -l accepts an integer (equal bins) or comma-separated weights
        if isinstance(machines, list):
            load_str = ",".join(str(x) for x in machines)
        else:
            load_str = str(machines)

        cmd = [
            str(self.exec_path), "makespan",
            "-i", str(input_file),
            "-o", str(output_file),
            "-k", str(field_cluster),
            "-m", str(method),
            "-w", str(weight_method),
            "-l", load_str,
        ]

        if class_columns is not None:
            cmd += ["-c", str(class_columns)]

        if extra_args:
            cmd += extra_args.split()

        return cmd

    def run_makespan(
        self,
        *,
        stats_file: Path | str,
        **build_kwargs: Any,
    ) -> Path:
        """Run ccphylo makespan, capturing stats from stdout.

        ccphylo makespan writes statistics to stdout. This method
        captures stdout and writes it to stats_file.

        Args:
            stats_file: Path to write statistics output.
            **build_kwargs: Arguments forwarded to build_command(),
                must include stats_file for build but it's handled here.

        Returns:
            Path to the stats file.
        """
        import subprocess

        build_kwargs["stats_file"] = stats_file
        cmd = self.build_command(**build_kwargs)
        stats_path = Path(stats_file)
        stats_path.parent.mkdir(parents=True, exist_ok=True)

        logger.info("Running: %s > %s", " ".join(cmd), stats_path)

        with open(stats_path, "w") as f:
            subprocess.run(
                cmd, stdout=f, stderr=subprocess.PIPE, check=True
            )

        return stats_path

    def map_outputs(
        self,
        workdir: Path,
        *,
        out_prefix: Optional[str] = None,
        app_args: Optional[Dict[str, Any]] = None,
        **kwargs: Any,
    ) -> Dict[str, Path]:
        if out_prefix is None:
            return {}
        return {
            "makespan": Path(str(out_prefix) + "_makespan.tsv"),
            "stats": Path(str(out_prefix) + "_makespan_stats.tsv"),
        }


class CCPhyloTreeApp(ApplicationRunner):
    """Runner for `ccphylo tree` — build a Newick tree from a distance matrix."""

    def __init__(self, exec_path: Path | str = "ccphylo") -> None:
        super().__init__(exec_path=Path(exec_path), tool_name="ccphylo")

    def build_command(
        self,
        *,
        input_file: Path | str,
        output_file: Path | str,
        method: str = "dnj",
        threads: int = 1,
        memory_disk: bool = False,
        tmp_dir: Optional[Path | str] = None,
        **kwargs: Any,
    ) -> List[str]:
        """Build `ccphylo tree` command.

        Args:
            input_file: Distance matrix (.phy).
            output_file: Output Newick tree file.
            method: Tree construction method (-m): "dnj", "nj", "upgma", etc.
            threads: Number of threads (-t).
            memory_disk: Allocate distance matrix on disk (-H).
            tmp_dir: Temporary directory (-T).
        """
        cmd = [
            str(self.exec_path), "tree",
            "-i", str(input_file),
            "-o", str(output_file),
            "-m", method,
            "-t", str(threads),
        ]

        if memory_disk:
            cmd.append("-H")

        if tmp_dir is not None:
            cmd += ["-T", str(tmp_dir)]

        return cmd

    def map_outputs(
        self,
        workdir: Path,
        *,
        out_prefix: Optional[str] = None,
        app_args: Optional[Dict[str, Any]] = None,
        **kwargs: Any,
    ) -> Dict[str, Path]:
        if out_prefix is None:
            return {}
        return {"tree": Path(str(out_prefix) + ".nwk")}
