"""MMseqs2 application runners for fast sequence searching and clustering."""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, List, Any
from cgecore.applications.base import ApplicationRunner, AlignerRunner
from cgecore.signatures.sample import SampleSpec
from cgecore.signatures.database import DatabaseSpec

import subprocess
import logging

logger = logging.getLogger(__name__)


class MMseqs2SearchApp(AlignerRunner):
    """Runner for `mmseqs easy-search` — fast sequence-vs-sequence search.

    Wraps mmseqs2's easy-search subcommand, which handles database creation
    internally.  For building a full distance matrix from a multi-FASTA
    (all-vs-all), see the pipeline's _distance_mmseqs2() method.

    Extends AlignerRunner with:
    - Sequence type detection (nucleotide / amino acid)
    - Configurable sensitivity and E-value
    - BLAST-tabular output parsing for identity scores
    """

    def __init__(self, exec_path: Path | str = "mmseqs") -> None:
        super().__init__(exec_path=Path(exec_path), tool_name="mmseqs")

    def _assemble_command(
        self,
        *,
        query: SampleSpec,
        db: DatabaseSpec,
        out_prefix: str,
        seq_type: str = "nucleotides",
        cores: int = 1,
        sensitivity: float = 7.5,
        evalue: float = 10000,
        min_seq_id: float = 0.0,
        tmp_dir: Optional[Path | str] = None,
        format_output: str = "query,target,fident,alnlen,mismatch,gapopen,"
                             "qstart,qend,tstart,tend,evalue,bits",
        **kwargs: Any,
    ) -> List[str]:
        """Build mmseqs easy-search command.

        Args:
            query: Query FASTA wrapped in SampleSpec.
            db: Target FASTA wrapped in DatabaseSpec (indexed=False).
                db.prefix is used as the target file path.
            out_prefix: Output file prefix. Result written to {out_prefix}.m8.
            seq_type: "nucleotides" (--search-type 3) or
                      "aminoacids" (--search-type 1).
            cores: Number of threads (--threads).
            sensitivity: Search sensitivity (-s), 1.0-7.5.
            evalue: E-value threshold (-e).
            min_seq_id: Minimum sequence identity filter (--min-seq-id).
            tmp_dir: Temporary directory for mmseqs.
            format_output: Output columns (--format-output).
        """
        if seq_type == "nucleotides":
            search_type = "3"
        elif seq_type == "aminoacids":
            search_type = "1"
        else:
            raise ValueError(
                f"seq_type must be 'nucleotides' or 'aminoacids', got '{seq_type}'"
            )

        output_file = str(out_prefix) + ".m8"
        tmp = str(tmp_dir) if tmp_dir is not None else str(Path(out_prefix).parent / "mmseqs_tmp")

        cmd = [
            str(self.exec_path), "easy-search",
            str(query.files[0]),
            str(db.prefix),
            output_file,
            tmp,
            "--search-type", search_type,
            "--threads", str(cores),
            "-s", str(sensitivity),
            "-e", str(evalue),
            "--format-output", format_output,
        ]

        if min_seq_id > 0:
            cmd += ["--min-seq-id", str(min_seq_id)]

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
        return {"result": Path(str(out_prefix) + ".m8")}

    def run_search(
        self,
        *,
        query_file: Path | str,
        target_file: Path | str,
        out_prefix: Path | str,
        seq_type: str = "nucleotides",
        cores: int = 1,
        sensitivity: float = 7.5,
        evalue: float = 10000,
        min_seq_id: float = 0.0,
        tmp_dir: Optional[Path | str] = None,
    ) -> Path:
        """Run mmseqs easy-search with raw file paths.

        Convenience method that accepts raw file paths instead of
        SampleSpec/DatabaseSpec.

        Args:
            query_file: Query FASTA file.
            target_file: Target FASTA file.
            out_prefix: Output prefix. Result written to {out_prefix}.m8.
            seq_type: "nucleotides" or "aminoacids".
            cores: Number of threads.
            sensitivity: Search sensitivity (1.0-7.5).
            evalue: E-value threshold.
            min_seq_id: Minimum sequence identity filter.
            tmp_dir: Temporary directory for mmseqs.

        Returns:
            Path to the .m8 result file.
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

        result = super(AlignerRunner, self).run(
            cmd=self.build_command(
                query=query_spec,
                db=db_spec,
                out_prefix=str(out_prefix),
                seq_type=seq_type,
                cores=cores,
                sensitivity=sensitivity,
                evalue=evalue,
                min_seq_id=min_seq_id,
                tmp_dir=tmp_dir,
            ),
            workdir=Path(out_prefix).parent,
            out_prefix=str(out_prefix),
        )

        return Path(str(out_prefix) + ".m8")

    @staticmethod
    def parse_result_file(
        result_file: Path | str,
        identity_col: int = 2,
    ) -> Dict[str, Dict[str, float]]:
        """Parse an mmseqs2 BLAST-tabular (.m8) output for identity scores.

        By default expects the standard BLAST-tab columns:
            query, target, fident, alnlen, mismatch, gapopen,
            qstart, qend, tstart, tend, evalue, bits

        For each query-target pair, keeps the hit with the highest identity.

        Args:
            result_file: Path to .m8 output file.
            identity_col: Column index (0-based) for fractional identity.
                Default 2 (fident column in standard BLAST-tab).

        Returns:
            Nested dict: {query: {target: {"identity": float, ...}}}.
        """
        results: Dict[str, Dict[str, Dict[str, float]]] = {}
        result_path = Path(result_file)

        if not result_path.is_file():
            return {}

        with open(result_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) <= identity_col:
                    continue

                query = parts[0]
                target = parts[1]

                try:
                    fident = float(parts[identity_col])
                except (ValueError, IndexError):
                    continue

                if query not in results:
                    results[query] = {}

                if target not in results[query] or fident > results[query][target]["identity"]:
                    results[query][target] = {"identity": fident}

        return results


class MMseqs2ClusterApp(ApplicationRunner):
    """Runner for `mmseqs easy-cluster` — fast sequence clustering.

    Wraps mmseqs2's easy-cluster subcommand, which clusters sequences
    by sequence identity without building a full distance matrix.
    """

    def __init__(self, exec_path: Path | str = "mmseqs") -> None:
        super().__init__(exec_path=Path(exec_path), tool_name="mmseqs")

    def build_command(
        self,
        *,
        input_file: Path | str,
        out_prefix: Path | str,
        seq_type: str = "nucleotides",
        min_seq_id: float = 0.3,
        threads: int = 1,
        tmp_dir: Optional[Path | str] = None,
        **kwargs: Any,
    ) -> List[str]:
        """Build mmseqs easy-cluster command.

        Args:
            input_file: Input FASTA file.
            out_prefix: Output prefix. Produces {out_prefix}_cluster.tsv.
            seq_type: "nucleotides" or "aminoacids".
            min_seq_id: Minimum sequence identity for clustering (0-1).
            threads: Number of threads.
            tmp_dir: Temporary directory for mmseqs.
        """
        if seq_type == "nucleotides":
            search_type = "3"
        elif seq_type == "aminoacids":
            search_type = "1"
        else:
            raise ValueError(
                f"seq_type must be 'nucleotides' or 'aminoacids', got '{seq_type}'"
            )

        tmp = str(tmp_dir) if tmp_dir is not None else str(
            Path(out_prefix).parent / "mmseqs_tmp"
        )

        cmd = [
            str(self.exec_path), "easy-cluster",
            str(input_file),
            str(out_prefix),
            tmp,
            "--search-type", search_type,
            "--min-seq-id", str(min_seq_id),
            "--threads", str(threads),
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
        return {"cluster_tsv": Path(str(out_prefix) + "_cluster.tsv")}

    def run_cluster(
        self,
        *,
        input_file: Path | str,
        out_prefix: Path | str,
        seq_type: str = "nucleotides",
        min_seq_id: float = 0.3,
        threads: int = 1,
        tmp_dir: Optional[Path | str] = None,
    ) -> Path:
        """Run mmseqs easy-cluster and return the cluster TSV path.

        Args:
            input_file: Input FASTA file.
            out_prefix: Output prefix.
            seq_type: "nucleotides" or "aminoacids".
            min_seq_id: Minimum sequence identity for clustering.
            threads: Number of threads.
            tmp_dir: Temporary directory for mmseqs.

        Returns:
            Path to the cluster TSV file ({out_prefix}_cluster.tsv).
        """
        cmd = self.build_command(
            input_file=input_file,
            out_prefix=out_prefix,
            seq_type=seq_type,
            min_seq_id=min_seq_id,
            threads=threads,
            tmp_dir=tmp_dir,
        )

        self.run(cmd=cmd, workdir=Path(out_prefix).parent)

        return Path(str(out_prefix) + "_cluster.tsv")

    @staticmethod
    def parse_cluster_tsv(path: Path | str) -> Dict[str, List[str]]:
        """Parse mmseqs2 easy-cluster TSV output.

        The file has two columns: representative_id and member_id.

        Returns:
            Dict mapping representative → list of member IDs.
        """
        clusters: Dict[str, List[str]] = defaultdict(list)
        cluster_path = Path(path)

        if not cluster_path.is_file():
            return {}

        with open(cluster_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                rep, member = parts[0], parts[1]
                clusters[rep].append(member)

        return dict(clusters)
