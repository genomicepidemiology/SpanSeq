"""SpanSeq command-line interface."""
from __future__ import annotations
import argparse
import logging
import os
import sys

from spanseq import __version__
from spanseq.config import SpanSeqConfig
from spanseq.pipeline import SpanSeqPipeline


def build_parser() -> argparse.ArgumentParser:
    """Build the SpanSeq argument parser."""
    parser = argparse.ArgumentParser(
        prog="spanseq",
        description="SpanSeq — homology-aware partitioning of biological sequence databases",
    )
    parser.add_argument(
        "-v", "--version", action="version", version=__version__,
    )

    subparsers = parser.add_subparsers(title="commands", dest="action")

    # ── Shared parent parser ──────────────────────────────────────
    parent = argparse.ArgumentParser(add_help=False)

    # Input
    input_grp = parent.add_argument_group("Input")
    input_mx = input_grp.add_mutually_exclusive_group(required=True)
    input_mx.add_argument(
        "-i", "--input_fasta", metavar="FILE",
        help="Multi-FASTA file with sequences",
    )
    input_mx.add_argument(
        "-if", "--input_folder", metavar="DIR",
        help="Folder with FASTA files",
    )
    input_mx.add_argument(
        "-ib", "--input_batch", metavar="FILE",
        help="File listing paths to FASTA files",
    )
    input_grp.add_argument(
        "-s", "--seqtype", required=True,
        choices=["nucleotides", "aminoacids"],
        help="Type of biological sequence",
    )

    # Output
    output_grp = parent.add_argument_group("Output")
    output_grp.add_argument(
        "-o", "--output_folder", required=True,
        help="Output directory",
    )
    output_grp.add_argument(
        "-f", "--output_format", default="minimal",
        choices=["minimal", "merged_table", "fasta_files"],
        help="Output format (default: minimal)",
    )
    output_grp.add_argument(
        "-tmp", "--temp_files", default=None,
        help="Location for temporary files",
    )
    output_grp.add_argument(
        "-r", "--keep_tmp", action="store_true", default=False,
        help="Do not remove temporary files",
    )

    # K-mer options
    kmer_grp = parent.add_argument_group("K-mer options")
    kmer_grp.add_argument("-k", "--kmer_size", type=int, default=None,
                          help="K-mer size for KMA or Mash")
    kmer_grp.add_argument("-m", "--minimizer_size", default=None,
                          help="Minimizer size for KMA")
    kmer_grp.add_argument("-p", "--prefix", default="-", type=str,
                          help="Sparse prefix for KMA (use 'TG' for large DBs)")
    kmer_grp.add_argument("-ME", "--MegaDB", action="store_true", default=False,
                          help="Enable KMA MegaDB mode for ~10^6 sequences")
    kmer_grp.add_argument("-l", "--max_length", type=int, default=None,
                          help="Approximate length of longest sequence (for Mash/GGSearch)")

    # Clustering
    clust_grp = parent.add_argument_group("Clustering / partitioning")
    clust_grp.add_argument("-c", "--min_dist", required=True, type=float,
                           help="Maximum distance for sequences in different bins (0-1)")
    clust_grp.add_argument("-b", "--bins", required=True,
                           help="Number of bins, or comma-separated proportions")

    # Makespan
    mspan_grp = parent.add_argument_group("Makespan options")
    mspan_grp.add_argument("-mP", "--makespanProcess", default="DBF",
                           choices=["DBF", "DFF"],
                           help="Makespan method (default: DBF)")
    mspan_grp.add_argument("-mW", "--makespanWeights", default="none",
                           choices=["none", "logX", "powX", "expX"],
                           help="Makespan weighting method")
    mspan_grp.add_argument("-mI", "--makespanImbalanced", default=None,
                           help="TSV file with class labels for imbalance-aware partitioning")

    # Tool paths
    tool_grp = parent.add_argument_group("Tool paths (optional, default: use PATH)")
    tool_grp.add_argument("-KP", "--kmaPath", default=None,
                          help="Path to KMA executable")
    tool_grp.add_argument("-CP", "--ccphyloPath", default=None,
                          help="Path to CCPhylo executable")
    tool_grp.add_argument("-MP", "--mashPath", default=None,
                          help="Path to Mash executable")
    tool_grp.add_argument("-DP", "--CDHitPath", default=None,
                          help="Path to CD-HIT executable")
    tool_grp.add_argument("-GP", "--GGSearchPath", default=None,
                          help="Path to ggsearch36 executable")
    tool_grp.add_argument("-SP", "--mmseqsPath", default=None,
                          help="Path to MMseqs2 executable")

    # Technical
    tech_grp = parent.add_argument_group("Technical")
    tech_grp.add_argument("-n", "--threads", default=1, type=int,
                          help="Number of threads")

    # ── Split subcommand ──────────────────────────────────────────
    split_parser = subparsers.add_parser(
        "split", parents=[parent],
        help="Partition sequences into bins with minimal cross-bin similarity",
    )
    split_grp = split_parser.add_argument_group("Split-specific options")
    split_grp.add_argument(
        "-d", "--distanceMethod", default="cosine",
        choices=["jaccard", "szymkiewicz_simpson", "cosine", "kmer_inv", "mash", "identity", "mmseqs2", "mmseqs-fast"],
        help="Distance calculation method (default: cosine)",
    )
    split_grp.add_argument(
        "-a", "--approach", default="all",
        choices=["all", "hobohm_reduce", "hobohm_split"],
        help="Pipeline mode (default: all)",
    )
    split_grp.add_argument(
        "-hd", "--hobohm1_distance", default=None, type=float,
        help="Hobohm1 identity threshold for pre-reduction",
    )
    split_grp.add_argument(
        "-hm", "--hobohm1_method", default="cdhit",
        choices=["kma", "cdhit"],
        help="Software for Hobohm1 reduction (default: cdhit)",
    )
    split_grp.add_argument(
        "-H", "--memory_disk", action="store_true",
        help="Allocate distance matrix on disk",
    )

    # ── Reduce subcommand ─────────────────────────────────────────
    subparsers.add_parser(
        "reduce", parents=[parent],
        help="Remove redundant sequences above a similarity threshold",
    )

    return parser


def main():
    """SpanSeq entry point."""
    parser = build_parser()
    args = parser.parse_args()

    if args.action is None:
        parser.print_help()
        sys.exit(1)

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    try:
        config = SpanSeqConfig.from_args(args)
        pipeline = SpanSeqPipeline(config)
        outputs = pipeline.run()

        print("\nSpanSeq completed successfully.")
        print("Output files:")
        for key, path in outputs.items():
            if isinstance(path, list):
                for p in path:
                    print(f"  {key}: {p}")
            else:
                print(f"  {key}: {path}")

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Configuration error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        logging.exception("Pipeline failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
