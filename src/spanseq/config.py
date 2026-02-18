"""SpanSeq configuration — maps CLI arguments to pipeline parameters."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Union
import math


# Distance methods available per tool
KMA_DISTANCE_METHODS = {
    "jaccard": 64,
    "szymkiewicz_simpson": 2048,
    "cosine": 256,
    "kmer_inv": 32,
}

KMA_DISTANCE_FACTORS = {
    "jaccard": 1.0,
    "szymkiewicz_simpson": 1.0,
    "cosine": 1.0,
    "kmer_inv": 100.0,
}

DISTANCE_TOOLS = {
    **{k: "kma" for k in KMA_DISTANCE_METHODS},
    "mash": "mash",
    "identity": "ggsearch36",
    "mmseqs2": "mmseqs2",
    "mmseqs-fast": "mmseqs-fast",
}


def compute_sketch_size(max_length: int) -> int:
    """Compute the smallest power of 2 >= max_length."""
    if max_length <= 0:
        raise ValueError("max_length must be positive")
    return 1 << math.ceil(math.log2(max_length))


@dataclass
class SpanSeqConfig:
    """Complete configuration for a SpanSeq run.

    Constructed from CLI arguments. All paths are resolved to absolute.
    """

    # Action
    action: str  # "split" or "reduce"

    # Input
    input_path: Path
    input_format: str  # "file", "folder", "batch"
    seq_type: str      # "nucleotides" or "aminoacids"

    # Output
    output_dir: Path
    output_format: str = "minimal"  # "minimal", "merged_table", "fasta_files"
    keep_tmp: bool = False

    # Distance / clustering
    min_dist: float = 0.3
    distance_method: str = "cosine"  # only for split
    approach: str = "all"            # "all", "hobohm_reduce", "hobohm_split"
    hobohm1_distance: Optional[float] = None
    hobohm1_method: str = "cdhit"    # "kma" or "cdhit"

    # Partitioning
    bins: Union[int, List[int]] = 5
    makespan_method: str = "DBF"       # "DBF" or "DFF"
    makespan_weights: str = "none"     # "none", "logX", "powX", "expX"
    imbalance_file: Optional[Path] = None
    class_columns: Optional[str] = None

    # K-mer / sketch
    kmer_size: Optional[int] = None
    minimizer_size: Optional[int] = None
    prefix: str = "-"
    megadb: bool = False
    max_length: Optional[int] = None
    sketch_size: Optional[int] = None

    # Tool paths (None = use PATH)
    kma_path: Optional[Path] = None
    mash_path: Optional[Path] = None
    cdhit_est_path: Optional[Path] = None
    cdhit_aa_path: Optional[Path] = None
    ccphylo_path: Optional[Path] = None
    ggsearch_path: Optional[Path] = None
    mmseqs_path: Optional[Path] = None

    # Technical
    threads: int = 1
    memory_disk: bool = False
    tmp_dir: Optional[Path] = None

    def __post_init__(self):
        """Resolve paths and compute derived values."""
        self.input_path = Path(self.input_path).resolve()
        self.output_dir = Path(self.output_dir).resolve()

        if self.tmp_dir is None:
            self.tmp_dir = self.output_dir / "tmp"
        else:
            self.tmp_dir = Path(self.tmp_dir).resolve()

        if self.imbalance_file is not None:
            self.imbalance_file = Path(self.imbalance_file).resolve()

        # Compute sketch size for mash
        if self.distance_method == "mash" and self.max_length is not None:
            self.sketch_size = compute_sketch_size(self.max_length)

    @property
    def distance_tool(self) -> str:
        """Which tool handles the selected distance method."""
        return DISTANCE_TOOLS[self.distance_method]

    @property
    def kma_dist_flag(self) -> Optional[int]:
        """KMA distance method flag, or None if not using KMA."""
        return KMA_DISTANCE_METHODS.get(self.distance_method)

    @property
    def kma_dist_factor(self) -> float:
        """Factor to divide min_dist by for KMA methods."""
        return KMA_DISTANCE_FACTORS.get(self.distance_method, 1.0)

    @property
    def effective_dist_value(self) -> float:
        """min_dist adjusted by the distance method's factor."""
        return self.min_dist / self.kma_dist_factor

    @property
    def results_dir(self) -> Path:
        return self.output_dir

    @property
    def log_dir(self) -> Path:
        return self.output_dir / "log"

    @property
    def sample_name(self) -> str:
        return self.input_path.stem

    @property
    def needs_hobohm(self) -> bool:
        return self.approach in ("hobohm_reduce", "hobohm_split")

    @classmethod
    def from_args(cls, args) -> SpanSeqConfig:
        """Build config from parsed argparse Namespace."""
        # Determine input format
        if getattr(args, "input_fasta", None) is not None:
            input_path = args.input_fasta
            input_format = "file"
        elif getattr(args, "input_folder", None) is not None:
            input_path = args.input_folder
            input_format = "folder"
        elif getattr(args, "input_batch", None) is not None:
            input_path = args.input_batch
            input_format = "batch"
        else:
            raise ValueError("No input specified (-i, -if, or -ib)")

        # Parse bins
        bins_raw = args.bins
        if "," in str(bins_raw):
            bins = [int(x) for x in str(bins_raw).split(",")]
        else:
            bins = int(bins_raw)

        kwargs = dict(
            action=args.action,
            input_path=input_path,
            input_format=input_format,
            seq_type=args.seqtype,
            output_dir=args.output_folder,
            output_format=getattr(args, "output_format", "minimal"),
            keep_tmp=getattr(args, "keep_tmp", False),
            min_dist=args.min_dist,
            bins=bins,
            makespan_method=getattr(args, "makespanProcess", "DBF"),
            makespan_weights=getattr(args, "makespanWeights", "none"),
            kmer_size=getattr(args, "kmer_size", None),
            minimizer_size=getattr(args, "minimizer_size", None) or None,
            prefix=getattr(args, "prefix", "-"),
            megadb=getattr(args, "MegaDB", False),
            max_length=getattr(args, "max_length", None),
            threads=int(getattr(args, "threads", 1)),
            memory_disk=getattr(args, "memory_disk", False),
            tmp_dir=getattr(args, "temp_files", None),
        )

        # Split-specific
        if args.action == "split":
            kwargs["distance_method"] = getattr(args, "distanceMethod", "cosine")
            kwargs["approach"] = getattr(args, "approach", "all")
            hd = getattr(args, "hobohm1_distance", None)
            kwargs["hobohm1_distance"] = float(hd) if hd else None
            kwargs["hobohm1_method"] = getattr(args, "hobohm1_method", "cdhit")

        # Imbalance file
        imb = getattr(args, "makespanImbalanced", None)
        if imb:
            kwargs["imbalance_file"] = imb

        # Tool paths
        for attr, key in [
            ("kmaPath", "kma_path"),
            ("mashPath", "mash_path"),
            ("CDHitPath", "cdhit_est_path"),
            ("ccphyloPath", "ccphylo_path"),
            ("GGSearchPath", "ggsearch_path"),
            ("mmseqsPath", "mmseqs_path"),
        ]:
            path = getattr(args, attr, None)
            if path:
                kwargs[key] = Path(path)

        return cls(**kwargs)
