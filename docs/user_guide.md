# SpanSeq User Guide

## Overview

SpanSeq is a tool for homology-aware partitioning of biological sequence databases. It splits a collection of DNA or protein sequences into bins (partitions) such that the maximum similarity between any two sequences in different bins stays below a user-defined threshold. This is essential for building unbiased training/validation/test sets in machine learning on biological sequences.

## Installation

### Prerequisites

SpanSeq depends on several external bioinformatics tools. The easiest way to install them is through conda:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Install from source

```bash
git clone https://github.com/genomicepidemiology/SpanSeq.git
cd SpanSeq

# Create environment with external tools
conda env create -n spanseq --file data/envs/spanseqenv.yml
conda activate spanseq

# Install SpanSeq
pip install -e .
```

### Verify installation

```bash
spanseq --version
spanseq split --help
```

## Quick Start

### Split a nucleotide database into 3 bins

```bash
spanseq split \
  -i sequences.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.7 \
  -b 3
```

This creates three files in `results/`:
- `sequences_clusters.tsv` — DBSCAN cluster assignments
- `sequences_makespan.tsv` — partition (bin) assignments
- `sequences_makespan_stats.tsv` — partition statistics

### Split with FASTA output files

```bash
spanseq split \
  -i sequences.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.7 \
  -b 3 \
  -f fasta_files
```

This additionally produces `sequences_M1.fsa`, `sequences_M2.fsa`, `sequences_M3.fsa` — one FASTA file per partition.

### Reduce redundancy

```bash
spanseq reduce \
  -i sequences.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.7 \
  -b 3
```

## Modes

### Split

Partitions a sequence database into *b* bins while ensuring that no two sequences in different bins share more than *c* similarity. The pipeline:

1. **(Optional)** Hobohm1 pre-reduction — removes highly similar sequences before distance calculation
2. **Distance matrix** — computes pairwise distances using the chosen method
3. **DBSCAN clustering** — groups sequences whose distance is below the threshold
4. **Makespan partitioning** — distributes clusters across bins to balance bin sizes

### Reduce

A simpler workflow that uses KMA's built-in Hobohm1 algorithm to cluster sequences above the similarity threshold, then partitions the clusters:

1. **KMA index + Hobohm1** — clusters sequences above the threshold
2. **Makespan partitioning** — distributes clusters across bins

## Distance Methods

SpanSeq supports several distance calculation methods. Choose the one that fits your data and accuracy needs:

| CLI flag (`-d`) | Tool | Sequence types | Description |
|---|---|---|---|
| `cosine` (default) | KMA | DNA | Cosine distance from k-mer profiles |
| `jaccard` | KMA | DNA | Jaccard distance from k-mer sets |
| `szymkiewicz_simpson` | KMA | DNA | Szymkiewicz-Simpson overlap coefficient |
| `kmer_inv` | KMA | DNA | Inverse k-mer frequency distance |
| `mash` | Mash | DNA, protein | MinHash-based Jaccard distance estimate |
| `identity` | GGSearch36 | DNA, protein | Global sequence identity (Smith-Waterman) |
| `mmseqs2` | MMseqs2 | DNA, protein | Sequence identity via MMseqs2 search |
| `mmseqs-fast` | MMseqs2 | DNA, protein | Direct clustering via MMseqs2 (skips distance matrix + DBSCAN) |

**When to use which:**

- **KMA methods** (`cosine`, `jaccard`, etc.) are fast and work well for large nucleotide databases. `cosine` is the default and recommended starting point.
- **Mash** is fast and memory-efficient for large datasets. Good for both DNA and protein. Requires `-l` (max sequence length).
- **Identity / GGSearch36** computes true global alignment identity. Most accurate but slowest (quadratic). Requires `-l`. Good for small-to-medium datasets.
- **MMseqs2** is a fast approximation of sequence identity. Scales well to large datasets. Good for both DNA and protein.
- **MMseqs2 fast** (`mmseqs-fast`) is the fastest option. It uses MMseqs2's native clustering directly, skipping the distance matrix and DBSCAN steps entirely. Best for very large datasets where speed is prioritized over fine-grained distance control.

## CLI Reference

### Common options (shared by split and reduce)

#### Input (mutually exclusive)
| Flag | Description |
|---|---|
| `-i`, `--input_fasta FILE` | Multi-FASTA file with sequences |
| `-if`, `--input_folder DIR` | Folder with FASTA files |
| `-ib`, `--input_batch FILE` | File listing paths to FASTA files |

#### Required
| Flag | Description |
|---|---|
| `-s`, `--seqtype {nucleotides,aminoacids}` | Type of biological sequence |
| `-o`, `--output_folder DIR` | Output directory |
| `-c`, `--min_dist FLOAT` | Maximum distance for sequences in different bins (0-1). Higher = more restrictive |
| `-b`, `--bins VALUE` | Number of bins (integer) or comma-separated proportions (e.g. `3,3,4`) |

#### Output options
| Flag | Default | Description |
|---|---|---|
| `-f`, `--output_format` | `minimal` | `minimal`: cluster/makespan TSVs only. `merged_table`: adds a merged partition table. `fasta_files`: adds per-partition FASTA files |
| `-tmp`, `--temp_files DIR` | `<output>/tmp` | Location for temporary files |
| `-r`, `--keep_tmp` | off | Keep temporary files after completion |

#### K-mer options
| Flag | Default | Description |
|---|---|---|
| `-k`, `--kmer_size INT` | auto | K-mer size for KMA or Mash |
| `-m`, `--minimizer_size INT` | auto | Minimizer size for KMA |
| `-p`, `--prefix STR` | `-` | Sparse prefix for KMA (use `TG` for large DBs) |
| `-ME`, `--MegaDB` | off | KMA MegaDB mode for ~10^6 sequences |
| `-l`, `--max_length INT` | — | Approximate length of longest sequence (required for Mash and GGSearch) |

#### Makespan options
| Flag | Default | Description |
|---|---|---|
| `-mP`, `--makespanProcess` | `DBF` | `DBF` (Decreasing Best Fit) or `DFF` (Decreasing First Fit) |
| `-mW`, `--makespanWeights` | `none` | Cluster weighting: `none`, `logX`, `powX`, `expX` |
| `-mI`, `--makespanImbalanced FILE` | — | TSV file with class labels for imbalance-aware partitioning |

#### Tool paths (optional — defaults to PATH lookup)
| Flag | Tool |
|---|---|
| `-KP`, `--kmaPath` | KMA |
| `-CP`, `--ccphyloPath` | CCPhylo |
| `-MP`, `--mashPath` | Mash |
| `-DP`, `--CDHitPath` | CD-HIT |
| `-GP`, `--GGSearchPath` | GGSearch36 |
| `-SP`, `--mmseqsPath` | MMseqs2 |

#### Technical
| Flag | Default | Description |
|---|---|---|
| `-n`, `--threads INT` | 1 | Number of threads |

### Split-specific options

| Flag | Default | Description |
|---|---|---|
| `-d`, `--distanceMethod` | `cosine` | Distance method (see table above) |
| `-a`, `--approach` | `all` | `all`: full all-vs-all. `hobohm_reduce`: pre-reduce, partition representatives only. `hobohm_split`: pre-reduce, then reassign clustered sequences |
| `-hd`, `--hobohm1_distance FLOAT` | — | Identity threshold for Hobohm1 pre-reduction (should be > `-c`) |
| `-hm`, `--hobohm1_method` | `cdhit` | `cdhit` or `kma` for Hobohm1 step |
| `-H`, `--memory_disk` | off | Allocate distance matrix on disk (for large datasets) |

## Output Files

### Minimal output (`-f minimal`, default)

| File | Description |
|---|---|
| `<sample>_clusters.tsv` | Tab-separated: sequence name, cluster ID, cluster weight. Produced by CCPhylo DBSCAN |
| `<sample>_makespan.tsv` | Tab-separated: sequence name, assigned partition. Produced by CCPhylo makespan |
| `<sample>_makespan_stats.tsv` | Partition statistics (sizes, balance metrics) |

### Merged table output (`-f merged_table`)

All of the above, plus:

| File | Description |
|---|---|
| `<sample>_partitions.tsv` | Merged table combining cluster and partition assignments per sequence |

### FASTA files output (`-f fasta_files`)

All of the above, plus:

| File | Description |
|---|---|
| `<sample>_M1.fsa` ... `<sample>_M<b>.fsa` | One FASTA file per partition containing the assigned sequences |

## Examples

### Protein sequences with Mash distance

```bash
spanseq split \
  -i proteins.fsa \
  -s aminoacids \
  -o results/ \
  -c 0.5 \
  -b 5 \
  -d mash \
  -l 1000 \
  -n 8
```

### Nucleotide sequences with Hobohm1 pre-reduction

```bash
spanseq split \
  -i genes.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.7 \
  -b 3 \
  -a hobohm_reduce \
  -hd 0.9 \
  -hm cdhit \
  -n 4
```

### Imbalance-aware partitioning

```bash
# labels.tsv: sequence_name<TAB>class_label
spanseq split \
  -i genes.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.7 \
  -b 3 \
  -mI labels.tsv
```

### MMseqs2 for large protein databases

```bash
spanseq split \
  -i large_proteins.fsa \
  -s aminoacids \
  -o results/ \
  -c 0.3 \
  -b 5 \
  -d mmseqs2 \
  -n 16
```

### Proportional bin sizes

```bash
# 60% training, 20% validation, 20% test
spanseq split \
  -i sequences.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.7 \
  -b 6,2,2
```

## References

If using SpanSeq, please cite:

> Ferrer Florensa, Alfred, et al. "SpanSeq: Similarity-based sequence data splitting method for improved development and assessment of deep learning projects." *NAR Genomics and Bioinformatics*, Volume 6, Issue 3, September 2024. https://doi.org/10.1093/nargab/lqae106

SpanSeq relies on these tools:

1. Clausen, Philip TLC, Frank M. Aarestrup, and Ole Lund. "Rapid and precise alignment of raw reads against redundant databases with KMA." *BMC bioinformatics* 19 (2018): 1-8.
2. Ondov, Brian D., et al. "Mash: fast genome and metagenome distance estimation using MinHash." *Genome biology* 17.1 (2016): 1-14.
3. Clausen, Philip TLC. "Scaling neighbor joining to one million taxa with dynamic and heuristic neighbor joining." *Bioinformatics* 39.1 (2023): btac774.
4. Hobohm, Uwe, et al. "Selection of representative protein data sets." *Protein Science* 1.3 (1992): 409-417.
