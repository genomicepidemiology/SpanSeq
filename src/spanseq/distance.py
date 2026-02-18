"""Distance computation methods and Hobohm1 reduction.

Provides a DistanceMixin class that SpanSeqPipeline inherits from.
Each method uses application runners (self.kma_index, self.mash, etc.)
and configuration (self.config) set up by the pipeline.
"""
from __future__ import annotations
from pathlib import Path
import logging
import subprocess

from spanseq.fasta import parse_fasta, write_fasta
from spanseq.applications.ggsearch import GGSearchApp
from spanseq.applications.mmseqs2 import MMseqs2SearchApp, MMseqs2ClusterApp

logger = logging.getLogger(__name__)


class DistanceMixin:
    """Mixin providing distance matrix computation and Hobohm1 reduction."""

    def _compute_distance(
        self, input_file: Path, output_file: Path
    ) -> Path:
        """Compute distance matrix using the configured tool."""
        cfg = self.config
        tool = cfg.distance_tool

        if tool == "kma":
            return self._distance_kma(input_file, output_file)
        elif tool == "mash":
            return self._distance_mash(input_file, output_file)
        elif tool == "ggsearch36":
            return self._distance_ggsearch(input_file, output_file)
        elif tool == "mmseqs2":
            return self._distance_mmseqs2(input_file, output_file)
        else:
            raise ValueError(f"Unknown distance tool: {tool}")

    def _distance_kma(self, input_file: Path, output_file: Path) -> Path:
        """KMA index + dist pipeline."""
        cfg = self.config
        index_prefix = cfg.tmp_dir / cfg.sample_name

        # Index
        index_cmd = self.kma_index.build_command(
            input_file=input_file,
            output_prefix=index_prefix,
            input_format=cfg.input_format,
            kmer_size=cfg.kmer_size,
            minimizer_size=cfg.minimizer_size,
            prefix=cfg.prefix,
            megadb=cfg.megadb,
            extra_args="",
        )
        self.kma_index.run(cmd=index_cmd, workdir=cfg.tmp_dir)

        # Distance
        dist_cmd = self.kma_dist.build_command(
            db_prefix=index_prefix,
            output_file=output_file,
            method=cfg.kma_dist_flag,
            threads=cfg.threads,
            tmp_dir=cfg.tmp_dir,
        )
        self.kma_dist.run(cmd=dist_cmd, workdir=cfg.tmp_dir)

        return output_file

    def _distance_mash(self, input_file: Path, output_file: Path) -> Path:
        """Mash triangle distance."""
        cfg = self.config

        if cfg.max_length is None:
            raise ValueError("--max_length is required when using mash distance")

        self.mash.run_triangle(
            output_file=output_file,
            input_file=input_file,
            input_format=cfg.input_format,
            seq_type=cfg.seq_type,
            kmer_size=cfg.kmer_size or 7,
            sketch_size=cfg.sketch_size or 1000,
            threads=cfg.threads,
        )
        return output_file

    def _distance_ggsearch(self, input_file: Path, output_file: Path) -> Path:
        """GGSearch36 all-vs-all distance matrix.

        Runs ggsearch36 incrementally: for each sequence, align it against
        all previously seen sequences, then append it to the target file.
        Builds a lower-triangular distance matrix in PHYLIP format.
        """
        cfg = self.config
        if cfg.max_length is None:
            raise ValueError("--max_length is required when using ggsearch distance")

        tmp = cfg.tmp_dir / "aln_files"
        tmp.mkdir(parents=True, exist_ok=True)

        target_file = tmp / "target.fsa"
        target_file.write_text("")  # empty target initially

        records = parse_fasta(input_file)
        n = len(records)
        seen_names = []

        with open(output_file, "w") as matrix_f:
            matrix_f.write(f"\t{n}\n")  # header: count

            for i, record in enumerate(records):
                name = record.id.replace("/", "-")
                query_file = tmp / f"{name}.fsa"

                # Write single query
                with open(query_file, "w") as qf:
                    write_fasta(qf, record)

                if i > 0:
                    # Align against target (all previous sequences)
                    result_file = tmp / f"{name}.ggsearch36"
                    self.ggsearch.run_alignment(
                        output_file=result_file,
                        query_file=query_file,
                        target_file=target_file,
                        seq_type=cfg.seq_type,
                        max_len=cfg.max_length,
                        cores=cfg.threads,
                    )
                    hits = GGSearchApp.parse_result_file(result_file)
                    query_hits = hits.get(record.id, {})

                    # Write matrix row
                    row = [name]
                    for prev_name in seen_names:
                        lookup = prev_name if len(prev_name) < 30 else prev_name[:31]
                        identity = query_hits.get(lookup, {}).get("identity", 0.0)
                        distance = 1.0 - identity
                        row.append(str(distance))
                    matrix_f.write("\t".join(row) + "\n")
                else:
                    matrix_f.write(f"{name}\n")

                # Append query to target for next iteration
                with open(target_file, "a") as tf:
                    write_fasta(tf, record)

                seen_names.append(name)
                logger.info("Aligned %d/%d sequences", i + 1, n)

        return output_file

    def _distance_mmseqs2(self, input_file: Path, output_file: Path) -> Path:
        """MMseqs2 all-vs-all identity distance matrix.

        Runs mmseqs easy-search of the input against itself, then converts
        the pairwise identity hits into a lower-triangular PHYLIP distance
        matrix (distance = 1 - identity).
        """
        cfg = self.config
        tmp = cfg.tmp_dir / "mmseqs_aln"
        tmp.mkdir(parents=True, exist_ok=True)

        out_prefix = tmp / cfg.sample_name
        self.mmseqs2.run_search(
            query_file=input_file,
            target_file=input_file,
            out_prefix=out_prefix,
            seq_type=cfg.seq_type,
            cores=cfg.threads,
            tmp_dir=tmp / "tmp",
        )

        # Parse results
        result_file = Path(str(out_prefix) + ".m8")
        hits = MMseqs2SearchApp.parse_result_file(result_file)

        # Get ordered sequence names from input
        records = parse_fasta(input_file)
        names = [r.id.replace("/", "-") for r in records]
        n = len(names)

        # Build lower-triangular distance matrix in PHYLIP format
        with open(output_file, "w") as f:
            f.write(f"\t{n}\n")
            for i, name_i in enumerate(names):
                if i == 0:
                    f.write(f"{name_i}\n")
                    continue
                row = [name_i]
                query_id = records[i].id
                query_hits = hits.get(query_id, {})
                for j in range(i):
                    target_id = records[j].id
                    hit = query_hits.get(target_id, {})
                    identity = hit.get("identity", 0.0)
                    distance = 1.0 - identity
                    row.append(str(distance))
                f.write("\t".join(row) + "\n")

        logger.info("MMseqs2 distance matrix written: %s", output_file)
        return output_file

    # ── Direct clustering (mmseqs-fast) ─────────────────────────────

    def _cluster_mmseqs2_fast(
        self, input_file: Path, clusters_file: Path
    ) -> Path:
        """Cluster sequences directly with MMseqs2 easy-cluster.

        Skips distance matrix and DBSCAN — uses MMseqs2's native clustering.
        Converts the output to DBSCAN-compatible format for makespan.

        Args:
            input_file: Input FASTA file.
            clusters_file: Path to write the clusters TSV.

        Returns:
            Path to clusters_file in DBSCAN-compatible format.
        """
        cfg = self.config
        tmp = cfg.tmp_dir / "mmseqs_cluster"
        tmp.mkdir(parents=True, exist_ok=True)

        out_prefix = tmp / cfg.sample_name
        min_seq_id = 1.0 - cfg.min_dist

        self.mmseqs2_cluster.run_cluster(
            input_file=input_file,
            out_prefix=out_prefix,
            seq_type=cfg.seq_type,
            min_seq_id=min_seq_id,
            threads=cfg.threads,
            tmp_dir=tmp / "tmp",
        )

        # Parse cluster output and convert to DBSCAN format
        cluster_tsv = Path(str(out_prefix) + "_cluster.tsv")
        raw_clusters = MMseqs2ClusterApp.parse_cluster_tsv(cluster_tsv)

        # Write DBSCAN-compatible format: #Sample\tNeighbors\tCluster
        with open(clusters_file, "w") as f:
            f.write("#Sample\tNeighbors\tCluster\n")
            for cluster_id, (rep, members) in enumerate(raw_clusters.items()):
                size = len(members)
                for member in members:
                    f.write(f"{member}\t{size}\t{cluster_id}\n")

        logger.info("MMseqs2 fast clustering: %d clusters written to %s",
                     len(raw_clusters), clusters_file)
        return clusters_file

    # ── Hobohm1 reduction ─────────────────────────────────────────────

    def _run_hobohm1(self, input_file: Path, output_prefix: Path) -> Path:
        """Run Hobohm1 redundancy reduction."""
        cfg = self.config

        if cfg.hobohm1_method == "cdhit":
            output_fasta = Path(str(output_prefix) + ".fsa")
            cmd = self.cdhit.build_command(
                input_file=input_file,
                output_file=output_fasta,
                seq_type=cfg.seq_type,
                threshold=1.0 - cfg.hobohm1_distance,
                threads=cfg.threads,
            )
            self.cdhit.run(cmd=cmd, workdir=cfg.tmp_dir)
            return output_fasta
        else:
            # KMA hobohm1 — run kma index with hobohm flags
            # Returns the original input (hobohm reduces via index)
            cmd = self.kma_index.build_command(
                input_file=input_file,
                output_prefix=output_prefix,
                input_format=cfg.input_format,
                kmer_size=cfg.kmer_size,
                minimizer_size=cfg.minimizer_size,
                hobohm1_value=1.0 - cfg.hobohm1_distance,
            )
            hobohm_clusters = Path(str(output_prefix) + "_hobohm.tsv")
            with open(hobohm_clusters, "w") as f:
                subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
            return input_file  # distance is still computed on original
