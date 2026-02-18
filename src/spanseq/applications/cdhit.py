"""CD-HIT application runner for sequence redundancy reduction."""
from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional, List, Any
from cgecore.applications.base import ApplicationRunner

import logging

logger = logging.getLogger(__name__)

# CD-HIT word size selection based on identity threshold.
# From the CD-HIT documentation.
NUCLEOTIDE_WORD_SIZES = [
    (0.90, 8),
    (0.88, 7),
    (0.85, 6),
    (0.80, 5),
    (0.75, 4),
]

PROTEIN_WORD_SIZES = [
    (0.70, 5),
    (0.60, 4),
    (0.50, 3),
    (0.40, 2),
]


def select_word_size(threshold: float, seq_type: str) -> int:
    """Select CD-HIT word size for a given identity threshold.

    Args:
        threshold: Sequence identity threshold (0-1).
        seq_type: "nucleotides" or "aminoacids".

    Returns:
        Appropriate word size.

    Raises:
        ValueError: If threshold is too low for any available word size.
    """
    table = NUCLEOTIDE_WORD_SIZES if seq_type == "nucleotides" else PROTEIN_WORD_SIZES
    for min_threshold, word_size in table:
        if threshold >= min_threshold:
            return word_size
    raise ValueError(
        f"Identity threshold {threshold} is too low for {seq_type}. "
        f"Minimum is {table[-1][0]}."
    )


class CdHitApp(ApplicationRunner):
    """Runner for CD-HIT / CD-HIT-EST sequence clustering.

    Automatically selects `cd-hit-est` for nucleotides and `cd-hit` for
    proteins, and picks the appropriate word size for the threshold.
    """

    def __init__(
        self,
        exec_path_est: Path | str = "cd-hit-est",
        exec_path_aa: Path | str = "cd-hit",
    ) -> None:
        # Validate at least one executable. We store both paths and pick
        # the right one at command build time.
        self._path_est = Path(exec_path_est)
        self._path_aa = Path(exec_path_aa)
        # Initialize with est as default for validation
        super().__init__(exec_path=self._path_est, tool_name="cd-hit")

    def build_command(
        self,
        *,
        input_file: Path | str,
        output_file: Path | str,
        seq_type: str = "nucleotides",
        threshold: float = 0.9,
        word_size: Optional[int] = None,
        threads: int = 1,
        extra_args: str = "",
        **kwargs: Any,
    ) -> List[str]:
        """Build CD-HIT command.

        Args:
            input_file: Input FASTA file.
            output_file: Output reduced FASTA.
            seq_type: "nucleotides" (cd-hit-est) or "aminoacids" (cd-hit).
            threshold: Sequence identity threshold (-c).
            word_size: Word size (-n). Auto-selected if None.
            threads: Number of threads (-T).
            extra_args: Additional CLI arguments.
        """
        if seq_type == "nucleotides":
            executable = str(self._path_est)
        else:
            executable = str(self._path_aa)

        if word_size is None:
            word_size = select_word_size(threshold, seq_type)

        cmd = [
            executable,
            "-i", str(input_file),
            "-o", str(output_file),
            "-c", str(threshold),
            "-n", str(word_size),
            "-T", str(threads),
        ]

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
            "fasta": base,
            "clstr": Path(str(base) + ".clstr"),
        }
