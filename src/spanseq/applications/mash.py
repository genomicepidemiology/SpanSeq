"""Mash distance application runner."""
from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional, List, Any
from cgecore.applications.base import ApplicationRunner

import subprocess
import logging

logger = logging.getLogger(__name__)


class MashApp(ApplicationRunner):
    """Runner for `mash triangle` — computes pairwise distance matrix.

    Mash writes its distance matrix to stdout, so this runner captures
    stdout and writes it to the output file.
    """

    def __init__(self, exec_path: Path | str = "mash") -> None:
        super().__init__(exec_path=Path(exec_path), tool_name="mash")

    def build_command(
        self,
        *,
        input_file: Path | str,
        input_format: str = "file",
        seq_type: str = "nucleotides",
        kmer_size: int = 7,
        sketch_size: int = 1000,
        threads: int = 1,
        extra_args_sketch: str = "",
        extra_args_triangle: str = "",
        **kwargs: Any,
    ) -> List[str]:
        """Build `mash triangle` command.

        Args:
            input_file: Input FASTA file or file list.
            input_format: "file" for -i, "batch" for -l.
            seq_type: "nucleotides" or "aminoacids". Aminoacids adds -a.
            kmer_size: K-mer size (-k).
            sketch_size: Sketch size (-s).
            threads: Number of threads (-p).
            extra_args_sketch: Extra sketch arguments.
            extra_args_triangle: Extra triangle arguments.
        """
        input_flag = "-l" if input_format == "batch" else "-i"

        cmd = [
            str(self.exec_path), "triangle",
            input_flag, str(input_file),
            "-k", str(kmer_size),
            "-s", str(sketch_size),
            "-i",
            "-p", str(threads),
        ]

        if seq_type == "aminoacids":
            cmd.append("-a")

        if extra_args_sketch:
            cmd += extra_args_sketch.split()

        if extra_args_triangle:
            cmd += extra_args_triangle.split()

        return cmd

    def run_triangle(
        self,
        *,
        output_file: Path | str,
        **build_kwargs: Any,
    ) -> Path:
        """Run mash triangle and write stdout to output_file.

        Args:
            output_file: Path where the distance matrix is written.
            **build_kwargs: Arguments forwarded to build_command().

        Returns:
            Path to the output distance matrix.

        Raises:
            subprocess.CalledProcessError: If mash exits with non-zero status.
        """
        cmd = self.build_command(**build_kwargs)
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        logger.info("Running: %s > %s", " ".join(cmd), output_path)

        with open(output_path, "w") as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)

        return output_path

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
