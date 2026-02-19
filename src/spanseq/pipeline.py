"""SpanSeq pipelines — split and reduce workflows as plain Python."""
from __future__ import annotations
from pathlib import Path
from typing import Optional
import logging
import shutil
import subprocess

from spanseq.config import SpanSeqConfig
from spanseq.distance import DistanceMixin
from spanseq.output import OutputMixin
from spanseq.applications.kma import KmaIndexApp, KmaDistApp
from spanseq.applications.mash import MashApp
from spanseq.applications.cdhit import CdHitApp
from spanseq.applications.ccphylo import CCPhyloDbscanApp, CCPhyloMakespanApp, CCPhyloTreeApp
from spanseq.applications.ggsearch import GGSearchApp
from spanseq.applications.mmseqs2 import MMseqs2SearchApp, MMseqs2ClusterApp

logger = logging.getLogger(__name__)


class SpanSeqPipeline(DistanceMixin, OutputMixin):
    """Orchestrates the SpanSeq split and reduce workflows.

    Pipeline steps for 'split':
        1. (Optional) Hobohm1 reduction via CD-HIT or KMA
        2. Distance matrix computation (KMA, Mash, or GGSearch36)
        3. DBSCAN clustering (CCPhylo)
        4. Makespan partitioning (CCPhylo)
        5. (Optional) Output formatting (merged table, FASTA files)

    Pipeline steps for 'reduce':
        1. KMA index with Hobohm1 → cluster output
        2. Makespan partitioning (CCPhylo)
    """

    def __init__(self, config: SpanSeqConfig) -> None:
        self.config = config
        self._setup_directories()
        self._init_apps()

    def _setup_directories(self):
        """Create output and temporary directories."""
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        self.config.tmp_dir.mkdir(parents=True, exist_ok=True)
        self.config.log_dir.mkdir(parents=True, exist_ok=True)

    def _resolve_tool(self, user_path: Optional[Path], default: str) -> Path:
        """Resolve a tool path: use user-provided path or find in PATH."""
        if user_path is not None:
            p = Path(user_path)
            if p.is_dir():
                # User gave a directory, look for the executable inside
                candidate = p / default
                if candidate.exists():
                    return candidate
                return p / default
            return p
        # Find in PATH
        found = shutil.which(default)
        if found is None:
            raise FileNotFoundError(
                f"'{default}' not found in PATH. "
                f"Install it or provide an explicit path."
            )
        return Path(found)

    def _init_apps(self):
        """Initialize application runners based on config."""
        cfg = self.config

        # Always need ccphylo for dbscan + makespan (+ tree if requested)
        ccphylo = self._resolve_tool(cfg.ccphylo_path, "ccphylo")
        self.dbscan = CCPhyloDbscanApp(exec_path=ccphylo)
        self.makespan = CCPhyloMakespanApp(exec_path=ccphylo)
        if cfg.tree:
            self.tree_app = CCPhyloTreeApp(exec_path=ccphylo)

        # Distance tool
        if cfg.action == "split":
            tool = cfg.distance_tool
            if tool == "kma":
                kma = self._resolve_tool(cfg.kma_path, "kma")
                self.kma_index = KmaIndexApp(exec_path=kma)
                self.kma_dist = KmaDistApp(exec_path=kma)
            elif tool == "mash":
                mash = self._resolve_tool(cfg.mash_path, "mash")
                self.mash = MashApp(exec_path=mash)
            elif tool == "ggsearch36":
                ggsearch = self._resolve_tool(cfg.ggsearch_path, "ggsearch36")
                self.ggsearch = GGSearchApp(exec_path=ggsearch)
            elif tool == "mmseqs2":
                mmseqs = self._resolve_tool(cfg.mmseqs_path, "mmseqs")
                self.mmseqs2 = MMseqs2SearchApp(exec_path=mmseqs)
            elif tool == "mmseqs-fast":
                mmseqs = self._resolve_tool(cfg.mmseqs_path, "mmseqs")
                self.mmseqs2_cluster = MMseqs2ClusterApp(exec_path=mmseqs)

            # Hobohm1 reduction tool (if needed)
            if cfg.needs_hobohm:
                if cfg.hobohm1_method == "cdhit":
                    est = self._resolve_tool(cfg.cdhit_est_path, "cd-hit-est")
                    aa = self._resolve_tool(cfg.cdhit_aa_path, "cd-hit")
                    self.cdhit = CdHitApp(exec_path_est=est, exec_path_aa=aa)
                else:
                    kma = self._resolve_tool(cfg.kma_path, "kma")
                    self.kma_index = KmaIndexApp(exec_path=kma)

        elif cfg.action == "reduce":
            kma = self._resolve_tool(cfg.kma_path, "kma")
            self.kma_index = KmaIndexApp(exec_path=kma)

    def run(self) -> dict:
        """Execute the pipeline based on config.action.

        Returns:
            Dict with paths to output files.
        """
        if self.config.action == "split":
            return self._run_split()
        elif self.config.action == "reduce":
            return self._run_reduce()
        else:
            raise ValueError(f"Unknown action: {self.config.action}")

    # ── Split pipeline ────────────────────────────────────────────────

    def _run_split(self) -> dict:
        """Execute the split pipeline."""
        cfg = self.config
        sample = cfg.sample_name
        tmp = cfg.tmp_dir

        logger.info("=== SpanSeq Split Pipeline ===")
        logger.info("Sample: %s", sample)
        logger.info("Distance method: %s (%s)", cfg.distance_method, cfg.distance_tool)

        clusters_file = cfg.results_dir / f"{sample}_clusters.tsv"

        if cfg.distance_method == "mmseqs-fast":
            # Fast path: MMseqs2 easy-cluster directly (no distance matrix / DBSCAN)
            input_for_cluster = cfg.input_path
            if cfg.needs_hobohm:
                input_for_cluster = self._run_hobohm1(
                    input_file=cfg.input_path,
                    output_prefix=tmp / f"{sample}_hobohm",
                )
            self._cluster_mmseqs2_fast(
                input_file=input_for_cluster,
                clusters_file=clusters_file,
            )
        else:
            # Standard path: distance matrix → DBSCAN
            input_for_distance = cfg.input_path
            if cfg.needs_hobohm:
                input_for_distance = self._run_hobohm1(
                    input_file=cfg.input_path,
                    output_prefix=tmp / f"{sample}_hobohm",
                )

            dist_file = self._compute_distance(
                input_file=input_for_distance,
                output_file=tmp / f"{sample}.phy",
            )
            logger.info("Distance matrix: %s", dist_file)

            # (Optional) Build tree from distance matrix
            if cfg.tree:
                tree_file = cfg.results_dir / f"{sample}.nwk"
                tree_cmd = self.tree_app.build_command(
                    input_file=dist_file,
                    output_file=tree_file,
                    method=cfg.tree_method,
                    threads=cfg.threads,
                    memory_disk=cfg.memory_disk,
                    tmp_dir=cfg.tmp_dir,
                )
                self.tree_app.run(cmd=tree_cmd, workdir=cfg.results_dir)
                logger.info("Tree: %s", tree_file)

            dbscan_cmd = self.dbscan.build_command(
                input_file=dist_file,
                output_file=clusters_file,
                dist_value=cfg.effective_dist_value,
                memory_disk=cfg.memory_disk,
                tmp_dir=cfg.tmp_dir,
            )
            self.dbscan.run(cmd=dbscan_cmd, workdir=cfg.results_dir)

        logger.info("Clusters: %s", clusters_file)

        # Step 3b (optional): Add class columns for imbalance-aware makespan
        makespan_input = clusters_file
        if cfg.imbalance_file is not None:
            makespan_input = self._add_class_columns(
                clusters_file=clusters_file,
                output_file=cfg.results_dir / f"{sample}_clusters_class.tsv",
            )

        # Step 4: Makespan partitioning
        makespan_file = cfg.results_dir / f"{sample}_makespan.tsv"
        stats_file = cfg.results_dir / f"{sample}_makespan_stats.tsv"
        self.makespan.run_makespan(
            input_file=makespan_input,
            output_file=makespan_file,
            stats_file=stats_file,
            machines=cfg.bins,
            field_cluster=3,
            method=cfg.makespan_method,
            weight_method=cfg.makespan_weights,
            class_columns=cfg.class_columns,
        )
        logger.info("Makespan: %s", makespan_file)

        outputs = {
            "clusters": clusters_file,
            "makespan": makespan_file,
            "stats": stats_file,
        }

        if cfg.tree:
            outputs["tree"] = tree_file

        # Step 5 (optional): Output formatting
        if cfg.output_format in ("merged_table", "fasta_files"):
            partitions_file = self._merge_tables(
                clusters_file=clusters_file,
                makespan_file=makespan_file,
                output_file=cfg.results_dir / f"{sample}_partitions.tsv",
            )
            outputs["partitions"] = partitions_file

        if cfg.output_format == "fasta_files":
            fasta_files = self._create_fasta_partitions(
                partitions_file=partitions_file,
                fasta_file=cfg.input_path,
            )
            outputs["fasta_files"] = fasta_files

        # Cleanup
        if not cfg.keep_tmp:
            logger.info("Removing temporary files: %s", cfg.tmp_dir)
            shutil.rmtree(cfg.tmp_dir, ignore_errors=True)

        logger.info("=== Pipeline complete ===")
        return outputs

    # ── Reduce pipeline ───────────────────────────────────────────────

    def _run_reduce(self) -> dict:
        """Execute the reduce pipeline."""
        cfg = self.config
        sample = cfg.sample_name
        tmp = cfg.tmp_dir

        logger.info("=== SpanSeq Reduce Pipeline ===")
        logger.info("Sample: %s", sample)

        # Step 1: KMA index with Hobohm1
        clusters_file = cfg.results_dir / f"{sample}_clusters.tsv"
        index_prefix = tmp / sample

        cmd = self.kma_index.build_command(
            input_file=cfg.input_path,
            output_prefix=index_prefix,
            input_format=cfg.input_format,
            kmer_size=cfg.kmer_size,
            minimizer_size=cfg.minimizer_size,
            prefix=cfg.prefix,
            megadb=cfg.megadb,
            hobohm1_value=cfg.min_dist,
        )

        # KMA hobohm1 writes cluster output to stdout
        with open(clusters_file, "w") as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
        logger.info("Clusters: %s", clusters_file)

        # Step 2: Makespan partitioning
        makespan_file = cfg.results_dir / f"{sample}_makespan.tsv"
        stats_file = cfg.results_dir / f"{sample}_makespan_stats.tsv"
        self.makespan.run_makespan(
            input_file=clusters_file,
            output_file=makespan_file,
            stats_file=stats_file,
            machines=cfg.bins,
            field_cluster=2,  # reduce uses column 2
            method=cfg.makespan_method,
            weight_method=cfg.makespan_weights,
        )
        logger.info("Makespan: %s", makespan_file)

        outputs = {
            "clusters": clusters_file,
            "makespan": makespan_file,
            "stats": stats_file,
        }

        if not cfg.keep_tmp:
            shutil.rmtree(cfg.tmp_dir, ignore_errors=True)

        logger.info("=== Pipeline complete ===")
        return outputs

