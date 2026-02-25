"""Output formatting methods for SpanSeq results.

Provides an OutputMixin class that SpanSeqPipeline inherits from.
Handles merging cluster/makespan tables, creating per-partition FASTA
files, and adding class columns for imbalance-aware partitioning.
"""
from __future__ import annotations
from pathlib import Path
import logging

from spanseq.fasta import parse_fasta, write_fasta

logger = logging.getLogger(__name__)


class OutputMixin:
    """Mixin providing output formatting methods."""

    def _merge_tables(
        self, clusters_file: Path, makespan_file: Path, output_file: Path
    ) -> Path:
        """Merge cluster and makespan tables into a single partition table.

        clusters TSV:  #Sample, Neighbors, Cluster   (sequence → cluster)
        makespan TSV:  #Cluster, Cluster_size, Cluster_weight, Partition
        output TSV:    id, cluster, partition         (sequence → partition)
        """
        import pandas as pd
        import io

        def _read(path: Path) -> pd.DataFrame:
            """Read a ccphylo TSV, skipping ## metadata comment lines."""
            with open(path) as f:
                content = "".join(line for line in f if not line.startswith("##"))
            return pd.read_csv(io.StringIO(content), sep="\t")

        clusters = _read(clusters_file)   # columns: #Sample, Neighbors, Cluster
        makespan = _read(makespan_file)   # columns: #Cluster, Cluster_size, Cluster_weight, Partition

        # Join on cluster ID: clusters.Cluster == makespan.#Cluster
        merged = pd.merge(
            clusters[["#Sample", "Cluster"]],
            makespan[["#Cluster", "Partition"]],
            left_on="Cluster",
            right_on="#Cluster",
            how="left",
        )
        merged = merged.rename(columns={
            "#Sample": "id",
            "Cluster": "cluster",
            "Partition": "partition",
        })
        merged[["id", "cluster", "partition"]].to_csv(output_file, sep="\t", index=False)

        logger.info("Merged table: %s", output_file)
        return output_file

    def _create_fasta_partitions(
        self, partitions_file: Path, fasta_file: Path
    ) -> list:
        """Split FASTA file into per-partition files."""
        import pandas as pd

        cfg = self.config
        df = pd.read_csv(partitions_file, sep="\t")

        bins_count = cfg.bins if isinstance(cfg.bins, int) else len(cfg.bins)
        sample = cfg.sample_name
        ext = fasta_file.suffix

        # Open output files
        output_files = {}
        for i in range(1, bins_count + 1):
            out_path = cfg.results_dir / f"{sample}_M{i}{ext}"
            output_files[i] = open(out_path, "w")

        try:
            for record in parse_fasta(fasta_file):
                match = df.loc[df["id"] == record.id, "partition"]
                if match.empty:
                    match = df.loc[df["id"] == record.description, "partition"]
                if match.empty:
                    logger.warning("Sequence %s not found in partition table", record.id)
                    continue
                partition = int(match.values[0])
                write_fasta(output_files[partition], record)
        finally:
            for f in output_files.values():
                f.close()

        result = [cfg.results_dir / f"{sample}_M{i}{ext}" for i in range(1, bins_count + 1)]
        logger.info("Created %d partition FASTA files", len(result))
        return result

    def _add_class_columns(self, clusters_file: Path, output_file: Path) -> Path:
        """Add class columns from imbalance file to cluster table."""
        import pandas as pd
        import io

        cfg = self.config

        # Skip ## metadata comment lines produced by ccphylo
        with open(clusters_file) as f:
            content = "".join(line for line in f if not line.startswith("##"))
        clusters = pd.read_csv(io.StringIO(content), sep="\t")

        classes = pd.read_csv(cfg.imbalance_file, sep="\t", header=None)

        # Rename class columns
        col_names = {0: "Names"}
        for i in range(1, len(classes.columns)):
            col_names[i] = f"Class_{i}"
        classes.rename(columns=col_names, inplace=True)

        merged = pd.merge(clusters, classes, how="inner",
                          left_on="#Sample", right_on="Names")
        del merged["Names"]

        if len(merged) < len(clusters):
            raise ValueError(
                "The imbalance file contains fewer sequences than the cluster table"
            )

        merged.to_csv(output_file, sep="\t", index=False)

        return output_file
