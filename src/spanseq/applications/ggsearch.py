"""GGSearch36 application runner for global-global sequence alignment."""
from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional, List, Any
from cgecore.applications.base import AlignerRunner
from cgecore.signatures.sample import SampleSpec
from cgecore.signatures.database import DatabaseSpec

import subprocess
import logging

logger = logging.getLogger(__name__)


class GGSearchApp(AlignerRunner):
    """Runner for `ggsearch36` — global-global sequence alignment.

    Extends AlignerRunner with:
    - Sequence type flag selection (nucleotide / amino acid)
    - stdout-to-file capture (ggsearch36 writes results to stdout)
    - Result file parsing for identity scores

    The standard AlignerRunner interface (query=SampleSpec, db=DatabaseSpec)
    is available via run(). For convenience, run_alignment() accepts raw
    file paths.
    """

    def __init__(self, exec_path: Path | str = "ggsearch36") -> None:
        super().__init__(exec_path=Path(exec_path), tool_name="ggsearch36")

    def _assemble_command(
        self,
        *,
        query: SampleSpec,
        db: DatabaseSpec,
        out_prefix: str,
        seq_type: str = "nucleotides",
        max_len: int = 1000,
        cores: int = 1,
        evalue: float = 20000,
        **kwargs: Any,
    ) -> List[str]:
        """Build ggsearch36 command.

        Args:
            query: Query FASTA wrapped in SampleSpec.
            db: Target FASTA wrapped in DatabaseSpec (indexed=False).
                db.prefix is used as the target file path.
            out_prefix: Output prefix (not used by ggsearch36 directly).
            seq_type: "nucleotides" (-n) or "aminoacids" (-p).
            max_len: Maximum sequence length for alignment (-M 1-max_len).
            cores: Number of threads (-T).
            evalue: E-value threshold (-E).
        """
        if seq_type == "nucleotides":
            seq_flag = "-n"
        elif seq_type == "aminoacids":
            seq_flag = "-p"
        else:
            raise ValueError(
                f"seq_type must be 'nucleotides' or 'aminoacids', got '{seq_type}'"
            )

        cmd = [
            str(self.exec_path),
            str(query.files[0]),
            str(db.prefix),
            seq_flag,
            "-m", "9i",
            "-T", str(cores),
            "-E", str(evalue),
            "-M", f"1-{max_len}",
            "-d", "0",
        ]

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
        return {"alignment": Path(str(out_prefix) + ".ggsearch36")}

    def run_alignment(
        self,
        *,
        output_file: Path | str,
        query_file: Path | str,
        target_file: Path | str,
        seq_type: str = "nucleotides",
        max_len: int = 1000,
        cores: int = 1,
        evalue: float = 20000,
    ) -> Path:
        """Run ggsearch36 and write output to file.

        Convenience method that accepts raw file paths instead of
        SampleSpec/DatabaseSpec.  ggsearch36 writes alignment results
        to stdout, so this method captures stdout and writes it to
        output_file.

        Args:
            output_file: Path where alignment output is written.
            query_file: Query FASTA file.
            target_file: Target FASTA file.
            seq_type: "nucleotides" or "aminoacids".
            max_len: Maximum sequence length.
            cores: Number of threads.
            evalue: E-value threshold.

        Returns:
            Path to the output file.
        """
        query_path = Path(query_file)
        target_path = Path(target_file)

        query_spec = SampleSpec(type="assembled", files=[query_path])
        db_spec = DatabaseSpec(
            name=target_path.name,
            root=target_path.parent,
            indexed=False,
            indexer=None,
            sequencetype="genes",
        )

        cmd = self.build_command(
            query=query_spec,
            db=db_spec,
            out_prefix=str(output_file),
            seq_type=seq_type,
            max_len=max_len,
            cores=cores,
            evalue=evalue,
        )

        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w") as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)

        return output_path

    @staticmethod
    def parse_result_file(ggsearch_file: Path | str) -> Dict[str, Dict[str, float]]:
        """Parse a ggsearch36 output file for identity scores.

        Args:
            ggsearch_file: Path to ggsearch36 output.

        Returns:
            Dict mapping subject names to
            {"identity": float, "score": float, "length": float}.
        """
        results = {}
        ggsearch_path = Path(ggsearch_file)

        if not ggsearch_path.is_file():
            return results

        with open(ggsearch_path, "r") as f:
            for line in f:
                if line.startswith("The best scores are:"):
                    for score_line in f:
                        if score_line.startswith(">>>"):
                            break
                        parts = score_line.strip().split()
                        if not parts:
                            continue
                        subject = parts[0]
                        # Parse from right: pos0=len, pos2=identity, pos5=score
                        non_empty = [p for p in reversed(parts) if p]
                        if len(non_empty) < 6:
                            continue

                        try:
                            len_aln = float(non_empty[0])
                            id_val = float(non_empty[2])
                            score_val = float(non_empty[5])
                        except (ValueError, IndexError):
                            continue

                        if subject not in results or id_val > results[subject]["identity"]:
                            results[subject] = {
                                "identity": id_val,
                                "score": score_val,
                                "length": len_aln,
                            }

        return results
