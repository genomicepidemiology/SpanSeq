"""KMA index and distance application runners."""
from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional, List, Any
from cgecore.applications.base import ApplicationRunner

import logging

logger = logging.getLogger(__name__)


class KmaIndexApp(ApplicationRunner):
    """Runner for `kma index` — builds k-mer index from FASTA input.

    When hobohm1_value is set, also performs Hobohm1 reduction and writes
    a cluster file to stdout (captured in the run result).
    """

    def __init__(self, exec_path: Path | str = "kma") -> None:
        super().__init__(exec_path=Path(exec_path), tool_name="kma")

    def build_command(
        self,
        *,
        input_file: Path | str,
        output_prefix: Path | str,
        input_format: str = "file",
        kmer_size: Optional[int] = None,
        minimizer_size: Optional[int] = None,
        prefix: str = "-",
        megadb: bool = False,
        hobohm1_value: Optional[float] = None,
        extra_args: str = "",
        **kwargs: Any,
    ) -> List[str]:
        """Build `kma index` command.

        Args:
            input_file: Input FASTA file or batch list.
            output_prefix: Output index prefix (e.g. tmp/sample).
            input_format: "file" for -i, "batch" for -batch.
            kmer_size: K-mer size (-k).
            minimizer_size: Minimizer size (-m). None to disable.
            prefix: Prefix option (-Sparse prefix). Default "-".
            megadb: Enable -ME flag for large databases.
            hobohm1_value: If set, enable Hobohm1 reduction (-hq/-ht value).
            extra_args: Additional CLI arguments as a string.
        """
        input_flag = "-batch" if input_format == "batch" else "-i"

        cmd = [
            str(self.exec_path), "index",
            input_flag, str(input_file),
            "-o", str(output_prefix),
        ]

        if kmer_size is not None:
            cmd += ["-k", str(kmer_size)]

        if hobohm1_value is not None:
            cmd += ["-hq", str(hobohm1_value), "-ht", str(hobohm1_value), "-and"]

        if minimizer_size is not None:
            cmd += ["-m", str(minimizer_size)]

        cmd += ["-Sparse", prefix]

        if megadb:
            cmd.append("-ME")

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
        base = Path(out_prefix)
        return {
            "comp": Path(str(base) + ".comp.b"),
            "length": Path(str(base) + ".length.b"),
            "name": Path(str(base) + ".name"),
            "seq": Path(str(base) + ".seq.b"),
        }


class KmaDistApp(ApplicationRunner):
    """Runner for `kma dist` — computes pairwise distance matrix."""

    def __init__(self, exec_path: Path | str = "kma") -> None:
        super().__init__(exec_path=Path(exec_path), tool_name="kma")

    def build_command(
        self,
        *,
        db_prefix: Path | str,
        output_file: Path | str,
        method: int = 256,
        threads: int = 1,
        tmp_dir: Optional[Path | str] = None,
        extra_args: str = "",
        **kwargs: Any,
    ) -> List[str]:
        """Build `kma dist` command.

        Args:
            db_prefix: KMA index prefix (from KmaIndexApp).
            output_file: Output distance matrix (.phy).
            method: Distance method flag (-d). Default 256 (cosine).
            threads: Number of threads (-t).
            tmp_dir: Temporary directory (-tmp).
            extra_args: Additional CLI arguments.
        """
        cmd = [
            str(self.exec_path), "dist",
            "-t_db", str(db_prefix),
            "-o", str(output_file),
            "-d", str(method),
            "-t", str(threads),
        ]

        if tmp_dir is not None:
            cmd += ["-tmp", str(tmp_dir)]

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
        return {"phy": Path(str(out_prefix) + ".phy")}
