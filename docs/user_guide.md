# SpanSeq User Guide

## Table of Contents

- [What is SpanSeq?](#what-is-spanseq)
- [Why do I need this?](#why-do-i-need-this)
- [Key concepts](#key-concepts)
- [Installation](#installation)
- [Tutorial: your first run](#tutorial-your-first-run)
- [Modes: split vs reduce](#modes-split-vs-reduce)
- [Choosing a distance method](#choosing-a-distance-method)
- [Understanding the parameters](#understanding-the-parameters)
- [Understanding the output files](#understanding-the-output-files)
- [Common use cases](#common-use-cases)
- [CLI reference](#cli-reference)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [References](#references)

---

## What is SpanSeq?

SpanSeq is a tool for **homology-aware partitioning** of biological sequence databases. Given a FASTA file of DNA or protein sequences, it splits them into separate bins (partitions) such that no two sequences in *different* bins are too similar to each other.

In practical terms: if you set a 70% similarity threshold, SpanSeq guarantees that every pair of sequences landing in different bins shares less than 70% similarity.

## Why do I need this?

When you build machine learning models on biological sequences (e.g. predicting antimicrobial resistance from gene sequences), you need to split your data into training, validation, and test sets. If you split randomly, closely related sequences can end up in both your training and test sets. The model then appears to perform well, but it is really just "memorizing" similar sequences it has already seen. This is called **data leakage**.

SpanSeq prevents data leakage by ensuring that similar sequences always stay in the *same* bin. When you use one bin for training and another for testing, you can be confident that the model is evaluated on genuinely different sequences.

**Example of the problem:**

Imagine you have two nearly identical gene variants (99% identity). If one is in your training set and the other in your test set, a model that simply memorizes sequences will get a "correct" prediction on the test sequence — not because it learned the biology, but because it saw an almost-identical copy during training. SpanSeq places both variants in the same bin, preventing this.

## Key concepts

Before using SpanSeq, it helps to understand a few terms:

### Similarity and distance

SpanSeq works with **distances** between sequences. Distance is the complement of similarity: if two sequences share 80% identity, their distance is 0.2 (i.e. `1 - 0.8`).

The `-c` (min_dist) parameter is the **distance threshold**. A value of `0.3` means: "ensure that all cross-bin sequence pairs have at least 30% distance" — equivalently, "at most 70% similarity".

| `-c` value | Max cross-bin similarity | Strictness |
|---|---|---|
| 0.1 | 90% | Very strict — only very distant sequences can be in different bins |
| 0.3 | 70% | Moderate — commonly used for gene-level analyses |
| 0.5 | 50% | Relaxed — sequences must be less than half identical |
| 0.7 | 30% | Very relaxed — allows moderately related sequences across bins |

**Higher `-c` = stricter partitioning = fewer, larger clusters = less balanced bins.**

### Bins and partitions

A **bin** (or partition) is one of the output groups. If you request 3 bins (`-b 3`), SpanSeq produces three groups of sequences. Typically:

- Bin 1 = training set
- Bin 2 = validation set
- Bin 3 = test set

You can also request **proportional bins** with `-b 6,2,2`, meaning 60% training, 20% validation, 20% test. SpanSeq will try to match these proportions as closely as possible while respecting the similarity constraint.

### Clusters

Before partitioning, SpanSeq groups similar sequences into **clusters**. All sequences in a cluster are more similar to each other than the threshold allows across bins. Entire clusters are then assigned to bins, which is what guarantees the similarity constraint.

### The pipeline

The split workflow follows these steps:

```
Input FASTA
    │
    ▼
(Optional) Hobohm1 reduction    ← removes near-identical sequences to speed things up
    │
    ▼
Distance matrix computation     ← measures how similar every pair of sequences is
    │
    ▼
DBSCAN clustering               ← groups similar sequences into clusters
    │
    ▼
Makespan partitioning           ← distributes clusters across bins, balancing sizes
    │
    ▼
Output files                    ← cluster assignments, bin assignments, statistics
```

## Installation

### Step 1: Install conda (if you don't have it)

SpanSeq depends on several bioinformatics tools that are easiest to install through conda. If you don't have conda or mamba, install [Miniforge](https://github.com/conda-forge/miniforge):

```bash
# Download miniforge (Linux)
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

### Step 2: Set up conda channels

SpanSeq's dependencies come from the bioconda and conda-forge channels. Make sure they are configured:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### Step 3: Clone and install SpanSeq

```bash
git clone https://github.com/genomicepidemiology/SpanSeq.git
cd SpanSeq

# Create a conda environment with all external tools
conda env create -n spanseq --file data/envs/spanseqenv.yml

# Activate the environment
conda activate spanseq

# Install SpanSeq itself
pip install -e .
```

### Step 4: Verify everything works

```bash
# Check SpanSeq is installed
spanseq --version

# Check the help for each command
spanseq split --help
spanseq reduce --help

# Verify that the external tools are available
kma -v
ccphylo -v
mash --version
cd-hit -h 2>&1 | head -1
```

If any of these commands fail, the corresponding tool is not installed. See [Troubleshooting](#troubleshooting).

### Updating SpanSeq

To update to the latest version:

```bash
cd SpanSeq
git pull
pip install -e .
```

## Tutorial: your first run

This tutorial walks you through a complete example step by step.

### Prepare your data

SpanSeq takes a multi-FASTA file as input. This is a text file where each sequence has a header line starting with `>` followed by the sequence:

```
>gene_001 beta-lactamase TEM-1
ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCC...
>gene_002 beta-lactamase TEM-2
ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCC...
>gene_003 carbapenemase NDM-1
ATGAATTGCCCAATATTATGCACCCGGTCGCGAA...
```

Save this file as (for example) `sequences.fsa`.

### Run SpanSeq split

Let's split the sequences into 3 bins with a 70% identity threshold (distance = 0.3):

```bash
spanseq split \
  -i sequences.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.3 \
  -b 3
```

Here's what each flag means:

| Flag | Value | Meaning |
|---|---|---|
| `-i` | `sequences.fsa` | Your input FASTA file |
| `-s` | `nucleotides` | These are DNA sequences (use `aminoacids` for proteins) |
| `-o` | `results/` | Output directory (created automatically) |
| `-c` | `0.3` | Distance threshold — sequences in different bins will have at most 70% similarity |
| `-b` | `3` | Split into 3 bins |

### What happens during the run

SpanSeq will print progress messages like:

```
2024-09-15 10:30:01 [INFO] spanseq.pipeline: === SpanSeq Split Pipeline ===
2024-09-15 10:30:01 [INFO] spanseq.pipeline: Sample: sequences
2024-09-15 10:30:01 [INFO] spanseq.pipeline: Distance method: cosine (kma)
2024-09-15 10:30:05 [INFO] spanseq.distance: Distance matrix: results/tmp/sequences.phy
2024-09-15 10:30:06 [INFO] spanseq.pipeline: Clusters: results/sequences_clusters.tsv
2024-09-15 10:30:06 [INFO] spanseq.pipeline: Makespan: results/sequences_makespan.tsv
2024-09-15 10:30:06 [INFO] spanseq.pipeline: === Pipeline complete ===
```

### Check the output

After the run, the `results/` directory contains:

```
results/
├── sequences_clusters.tsv         ← which cluster each sequence belongs to
├── sequences_makespan.tsv         ← which bin each sequence is assigned to
└── sequences_makespan_stats.tsv   ← statistics about bin sizes
```

Open `sequences_makespan.tsv` to see the bin assignments:

```
#Sample    Machine
gene_001   1
gene_002   1
gene_003   2
gene_004   3
...
```

The `Machine` column is the bin number (1, 2, or 3). You can use this to create your train/validation/test split.

### Get FASTA output files directly

If you want SpanSeq to produce the actual FASTA files for each bin (instead of just the assignment table), add `-f fasta_files`:

```bash
spanseq split \
  -i sequences.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.3 \
  -b 3 \
  -f fasta_files
```

This produces everything above, plus:

```
results/
├── ...
├── sequences_M1.fsa   ← all sequences assigned to bin 1
├── sequences_M2.fsa   ← all sequences assigned to bin 2
└── sequences_M3.fsa   ← all sequences assigned to bin 3
```

You can load these directly in your ML pipeline.

## Modes: split vs reduce

SpanSeq has two modes: **split** and **reduce**.

### Split mode

**Use split when:** you want to partition your data for machine learning (training/validation/test sets).

Split performs a full analysis:
1. Computes pairwise distances between all sequences
2. Clusters sequences that are too similar to be in different bins
3. Distributes clusters across bins, balancing their sizes

```bash
spanseq split -i sequences.fsa -s nucleotides -o results/ -c 0.3 -b 3
```

### Reduce mode

**Use reduce when:** you want to remove redundancy from your dataset by clustering similar sequences, then partition the remaining representatives.

Reduce is simpler and faster:
1. Uses KMA's Hobohm1 algorithm to cluster sequences above the threshold
2. Distributes clusters across bins

```bash
spanseq reduce -i sequences.fsa -s nucleotides -o results/ -c 0.3 -b 3
```

**Key difference:** Split computes a full distance matrix and uses DBSCAN clustering, which is more accurate. Reduce uses a greedy algorithm (Hobohm1), which is faster but less precise.

### When to use which?

| Scenario | Recommended mode |
|---|---|
| Preparing ML train/val/test splits | `split` |
| Small-to-medium dataset (< 50k sequences) | `split` |
| Very large dataset (> 100k sequences), speed matters | `reduce` or `split -d mmseqs-fast` |
| Just want to remove redundancy | `reduce` |

## Choosing a distance method

The distance method (`-d` flag) determines how SpanSeq measures similarity between sequences. This is the most important choice after the threshold.

### Available methods

| Method | Flag | Tool | Works with | Speed | Accuracy |
|---|---|---|---|---|---|
| Cosine distance | `-d cosine` | KMA | DNA only | Fast | Good |
| Jaccard distance | `-d jaccard` | KMA | DNA only | Fast | Good |
| Szymkiewicz-Simpson | `-d szymkiewicz_simpson` | KMA | DNA only | Fast | Good |
| Inverse k-mer frequency | `-d kmer_inv` | KMA | DNA only | Fast | Good |
| Mash distance | `-d mash` | Mash | DNA + protein | Fast | Good |
| Global identity | `-d identity` | GGSearch36 | DNA + protein | Slow | Best |
| MMseqs2 identity | `-d mmseqs2` | MMseqs2 | DNA + protein | Medium | Good |
| MMseqs2 fast cluster | `-d mmseqs-fast` | MMseqs2 | DNA + protein | Fastest | Approximate |

### Decision guide

**Start here:**

1. **DNA sequences?**
   - For most cases, use the default (`-d cosine`). It is fast and works well.
   - If you need true sequence identity (e.g. for publication), use `-d identity` (but it is much slower).

2. **Protein sequences?**
   - Use `-d mmseqs2` for a good balance of speed and accuracy.
   - Use `-d mash` if you have very many sequences (> 50k).
   - Use `-d identity` for the most accurate results on small datasets (< 10k sequences).

3. **Very large dataset (> 100k sequences)?**
   - Use `-d mmseqs-fast` — it skips the distance matrix entirely and uses MMseqs2's built-in clustering. This is the fastest option by far, but gives you less control over the exact distance threshold.

4. **Need exact global alignment identity?**
   - Use `-d identity`. This runs GGSearch36, which performs Smith-Waterman global alignment between every pair of sequences. It is the gold standard for accuracy but has quadratic runtime — impractical for more than ~10k sequences.

### Important notes about specific methods

**KMA methods** (`cosine`, `jaccard`, `szymkiewicz_simpson`, `kmer_inv`):
- Work only with DNA sequences, not protein.
- Use k-mer-based distance approximations.
- The default k-mer size is determined automatically by KMA.
- `cosine` is the default and recommended for most nucleotide datasets.

**Mash** (`-d mash`):
- Requires the `-l` flag (approximate length of the longest sequence).
- Uses MinHash sketching to approximate Jaccard distance.
- Very memory-efficient — handles hundreds of thousands of sequences.

**GGSearch36** (`-d identity`):
- Requires the `-l` flag.
- Computes true global alignment identity between all pairs.
- Runtime is O(n^2) — scales poorly beyond ~10k sequences.
- Most biologically meaningful if you care about percent identity thresholds.

**MMseqs2** (`-d mmseqs2`):
- Fast approximate sequence identity using prefiltering + alignment.
- Good balance of speed and accuracy.
- Works well up to ~100k sequences.

**MMseqs2 fast** (`-d mmseqs-fast`):
- Uses `mmseqs easy-cluster` directly instead of computing a full distance matrix.
- Skips the DBSCAN step entirely.
- Fastest option available, suitable for hundreds of thousands of sequences.
- The clustering threshold is approximate — if you need precise control over the identity threshold, use `-d mmseqs2` instead.

## Understanding the parameters

### The distance threshold (`-c`)

This is the most important parameter. It controls how strict the partitioning is.

`-c 0.3` means: "the maximum allowed distance between sequences in different bins is 0.3." Since distance = 1 - similarity, this is equivalent to saying "sequences in different bins share at most 70% similarity."

**How to choose a good threshold:**

- There is no universal "correct" value — it depends on your data and your question.
- A common starting point for gene-level analyses is `-c 0.3` (70% identity).
- For closely related sequences (e.g. variants of the same gene), you may need a stricter threshold like `-c 0.1` (90% identity).
- For protein families with remote homology, a relaxed threshold like `-c 0.5` or higher may be appropriate.
- **Start with a moderate value and inspect the output.** If bins are very unbalanced (one bin has 90% of the sequences), try a less strict threshold. If the makespan stats show good balance, your threshold is reasonable.

### Number of bins (`-b`)

The simplest usage is an integer:

```bash
-b 3    # Three equally-sized bins
-b 5    # Five equally-sized bins
```

For unequal splits (e.g. 80/10/10 for train/val/test), use comma-separated proportions:

```bash
-b 8,1,1    # 80% / 10% / 10%
-b 6,2,2    # 60% / 20% / 20%
-b 7,3      # 70% / 30% (two bins)
```

The numbers are **relative weights**, not percentages. `-b 3,3,4` and `-b 30,30,40` produce the same result.

**Note:** bins may not be perfectly balanced because entire clusters must go into the same bin. The makespan algorithm does its best to balance sizes given this constraint.

### Hobohm1 pre-reduction (`-a`, `-hd`, `-hm`)

For large datasets, you can speed up the distance computation by first removing near-identical sequences. This is the Hobohm1 algorithm.

The `-a` flag controls the approach:

- `-a all` (default): no pre-reduction, compute distances between all sequences.
- `-a hobohm_reduce`: pre-reduce the dataset, compute distances only between representatives, and partition the representatives. Clustered sequences are not reassigned.
- `-a hobohm_split`: same as `hobohm_reduce`, but after partitioning the representatives, the clustered sequences are reassigned to the same bins as their representatives.

If you use Hobohm1, you also need:

- `-hd`: the identity threshold for the pre-reduction. This should be **stricter (higher) than `-c`**. For example, if `-c 0.3`, use `-hd 0.1` (meaning: collapse sequences with > 90% identity before the main analysis).
- `-hm`: which tool to use for the reduction (`cdhit` or `kma`). CD-HIT is the default and recommended.

**Example:**

```bash
spanseq split \
  -i large_dataset.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.3 \
  -b 3 \
  -a hobohm_split \
  -hd 0.1 \
  -hm cdhit \
  -n 8
```

This first collapses sequences sharing > 90% identity, then runs the full pipeline on the remaining representatives, then assigns all original sequences back to bins.

### Threading (`-n`)

Most distance methods support multi-threading. Set `-n` to the number of CPU cores you want to use:

```bash
-n 8    # Use 8 threads
```

The speedup depends on the distance method and dataset size. For KMA and MMseqs2, multi-threading provides significant acceleration. For GGSearch36, the speedup is more modest.

### Output format (`-f`)

| Value | What you get |
|---|---|
| `minimal` (default) | Cluster assignments + bin assignments + statistics (3 TSV files) |
| `merged_table` | Everything above + a merged table joining clusters and bins |
| `fasta_files` | Everything above + one FASTA file per bin |

Use `fasta_files` if you want ready-to-use FASTA files for each bin. Use `minimal` if you just need the assignment tables and will do the splitting yourself.

### Memory disk mode (`-H`)

For very large datasets, the distance matrix may not fit in RAM. The `-H` flag tells CCPhylo to allocate the distance matrix on disk instead of in memory:

```bash
spanseq split -i huge_dataset.fsa -s nucleotides -o results/ -c 0.3 -b 3 -H
```

This is slower but allows you to process datasets that would otherwise run out of memory.

## Understanding the output files

### Clusters file (`<sample>_clusters.tsv`)

Tab-separated file produced by DBSCAN clustering:

```
#Sample    Neighbors    Cluster
gene_001   5            0
gene_002   5            0
gene_003   5            0
gene_004   3            1
gene_005   3            1
gene_006   1            2
```

| Column | Meaning |
|---|---|
| `#Sample` | Sequence name (from the FASTA header) |
| `Neighbors` | Number of sequences in this cluster |
| `Cluster` | Cluster ID (sequences in the same cluster are similar) |

### Makespan file (`<sample>_makespan.tsv`)

Tab-separated file with the final bin assignment:

```
#Sample    Machine
gene_001   1
gene_002   1
gene_003   1
gene_004   2
gene_005   2
gene_006   3
```

| Column | Meaning |
|---|---|
| `#Sample` | Sequence name |
| `Machine` | Bin number (1-indexed) |

This is the file you use to create your training/validation/test splits.

### Statistics file (`<sample>_makespan_stats.tsv`)

Summary statistics about the partitioning. Shows how many sequences ended up in each bin and balance metrics. Use this to check if your bins are reasonably balanced.

### Partitions file (`<sample>_partitions.tsv`, with `-f merged_table`)

A merged table combining the cluster and makespan information:

```
#Sample    Neighbors    Cluster    Machine
gene_001   5            0          1
gene_002   5            0          1
...
```

### Per-bin FASTA files (`<sample>_M1.fsa`, etc., with `-f fasta_files`)

One FASTA file per bin, containing all sequences assigned to that bin. Ready to use in your ML pipeline.

## Common use cases

### 1. Standard ML train/val/test split (DNA)

The most common use case. Split nucleotide sequences into 3 bins at 70% identity:

```bash
spanseq split \
  -i genes.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.3 \
  -b 3 \
  -f fasta_files \
  -n 4
```

**Output:** `genes_M1.fsa`, `genes_M2.fsa`, `genes_M3.fsa`.
Use bin 1 for training, bin 2 for validation, bin 3 for testing.

### 2. Custom proportions (80/10/10)

```bash
spanseq split \
  -i genes.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.3 \
  -b 8,1,1 \
  -f fasta_files \
  -n 4
```

### 3. Protein sequences with MMseqs2

```bash
spanseq split \
  -i proteins.fsa \
  -s aminoacids \
  -o results/ \
  -c 0.3 \
  -b 3 \
  -d mmseqs2 \
  -f fasta_files \
  -n 8
```

### 4. Large dataset with Mash

For large datasets (> 50k sequences), Mash is fast and memory-efficient. You need to provide the approximate length of the longest sequence with `-l`:

```bash
spanseq split \
  -i large_dataset.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.3 \
  -b 3 \
  -d mash \
  -l 3000 \
  -n 8
```

### 5. Very large dataset with MMseqs2 fast clustering

For very large datasets (> 100k sequences) where speed is the priority:

```bash
spanseq split \
  -i huge_dataset.fsa \
  -s aminoacids \
  -o results/ \
  -c 0.3 \
  -b 5 \
  -d mmseqs-fast \
  -n 16
```

### 6. Accurate global identity with GGSearch36

For small datasets (< 10k sequences) where you need true global alignment identity:

```bash
spanseq split \
  -i small_dataset.fsa \
  -s aminoacids \
  -o results/ \
  -c 0.3 \
  -b 3 \
  -d identity \
  -l 500 \
  -n 4
```

### 7. Large dataset with Hobohm1 pre-reduction

When you have many highly similar sequences (e.g. from public databases), pre-reduce them before the main analysis:

```bash
spanseq split \
  -i big_dataset.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.3 \
  -b 3 \
  -a hobohm_split \
  -hd 0.1 \
  -hm cdhit \
  -n 8
```

This first collapses sequences above 90% identity, runs the full pipeline on the remaining representatives, and then assigns all original sequences to the same bin as their representative.

### 8. Imbalance-aware partitioning

When your sequences have class labels (e.g. "resistant" / "susceptible"), SpanSeq can try to balance the class distribution across bins:

First, create a labels file (`labels.tsv`):

```
gene_001	resistant
gene_002	resistant
gene_003	susceptible
gene_004	susceptible
```

Then run:

```bash
spanseq split \
  -i genes.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.3 \
  -b 3 \
  -mI labels.tsv
```

### 9. Redundancy reduction

Simply reduce redundancy without ML partitioning in mind:

```bash
spanseq reduce \
  -i sequences.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.3 \
  -b 3
```

### 10. Custom tool paths

If the external tools are not in your PATH, specify their locations explicitly:

```bash
spanseq split \
  -i sequences.fsa \
  -s nucleotides \
  -o results/ \
  -c 0.3 \
  -b 3 \
  -KP /opt/kma/kma \
  -CP /opt/ccphylo/ccphylo
```

## CLI reference

### Global options

```bash
spanseq --version        # Print version
spanseq --help           # Show help
spanseq split --help     # Show split-specific help
spanseq reduce --help    # Show reduce-specific help
```

### Common options (shared by split and reduce)

#### Input (mutually exclusive — exactly one is required)

| Flag | Description |
|---|---|
| `-i`, `--input_fasta FILE` | Multi-FASTA file with sequences |
| `-if`, `--input_folder DIR` | Folder containing FASTA files (processed as a single merged input) |
| `-ib`, `--input_batch FILE` | Text file listing paths to FASTA files (one per line) |

#### Required options

| Flag | Description |
|---|---|
| `-s`, `--seqtype` | Sequence type: `nucleotides` or `aminoacids` |
| `-o`, `--output_folder DIR` | Output directory (created if it doesn't exist) |
| `-c`, `--min_dist FLOAT` | Distance threshold (0 to 1). See [Understanding the parameters](#the-distance-threshold--c) |
| `-b`, `--bins VALUE` | Number of bins (integer) or proportions (comma-separated, e.g. `6,2,2`) |

#### Output options

| Flag | Default | Description |
|---|---|---|
| `-f`, `--output_format` | `minimal` | Output format: `minimal`, `merged_table`, or `fasta_files` |
| `-tmp`, `--temp_files DIR` | `<output>/tmp` | Location for temporary files |
| `-r`, `--keep_tmp` | off | Keep temporary files after completion (useful for debugging) |

#### K-mer options

| Flag | Default | Description |
|---|---|---|
| `-k`, `--kmer_size INT` | auto | K-mer size for KMA or Mash. Usually you don't need to set this |
| `-m`, `--minimizer_size INT` | auto | Minimizer size for KMA indexing |
| `-p`, `--prefix STR` | `-` | Sparse prefix for KMA (use `TG` for large databases) |
| `-ME`, `--MegaDB` | off | Enable KMA MegaDB mode for datasets with ~1 million sequences |
| `-l`, `--max_length INT` | — | Approximate length of the longest sequence. **Required** for Mash (`-d mash`) and GGSearch (`-d identity`) |

#### Makespan options

| Flag | Default | Description |
|---|---|---|
| `-mP`, `--makespanProcess` | `DBF` | Makespan scheduling algorithm: `DBF` (Decreasing Best Fit) or `DFF` (Decreasing First Fit) |
| `-mW`, `--makespanWeights` | `none` | Cluster weighting: `none`, `logX`, `powX`, or `expX` |
| `-mI`, `--makespanImbalanced FILE` | — | TSV file with class labels for imbalance-aware partitioning |

#### Tool paths (optional — by default SpanSeq finds tools in your PATH)

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
| `-n`, `--threads INT` | 1 | Number of CPU threads to use |

### Split-specific options

These options are only available with `spanseq split`:

| Flag | Default | Description |
|---|---|---|
| `-d`, `--distanceMethod` | `cosine` | Distance method. See [Choosing a distance method](#choosing-a-distance-method) |
| `-a`, `--approach` | `all` | `all`: full all-vs-all. `hobohm_reduce`: pre-reduce, partition representatives. `hobohm_split`: pre-reduce, then reassign |
| `-hd`, `--hobohm1_distance FLOAT` | — | Identity threshold for Hobohm1 pre-reduction. Should be stricter than `-c` |
| `-hm`, `--hobohm1_method` | `cdhit` | Tool for Hobohm1 reduction: `cdhit` or `kma` |
| `-H`, `--memory_disk` | off | Allocate distance matrix on disk instead of RAM (for very large datasets) |

## Troubleshooting

### "X not found in PATH"

```
Error: 'kma' not found in PATH. Install it or provide an explicit path.
```

This means SpanSeq can't find one of the external tools. Solutions:

1. **Activate your conda environment:** `conda activate spanseq`
2. **Check the tool is installed:** `which kma` (replace `kma` with the tool name)
3. **Provide the path explicitly:** `-KP /path/to/kma`
4. **Reinstall the tool:** `conda install -c bioconda kma`

### "--max_length is required"

```
ValueError: --max_length is required when using mash distance
```

Mash and GGSearch36 require the `-l` flag. Set it to the approximate length of your longest sequence:

```bash
# Check your longest sequence
grep -v ">" sequences.fsa | awk '{print length}' | sort -rn | head -1
# Use that value (rounded up) as -l
```

### Very unbalanced bins

If most sequences end up in one bin, your threshold (`-c`) is probably too strict for your dataset. Try:

1. **Lowering `-c`:** e.g. from 0.3 to 0.5
2. **Using more bins:** e.g. from 3 to 5
3. **Using Hobohm1 pre-reduction:** `-a hobohm_split -hd 0.1`

### Out of memory

For very large datasets, the distance matrix can exhaust your RAM. Options:

1. **Use disk mode:** add `-H` to store the matrix on disk
2. **Use a faster method:** `-d mmseqs-fast` skips the distance matrix entirely
3. **Pre-reduce your data:** use `-a hobohm_split` to collapse redundant sequences first
4. **Use Mash:** `-d mash` is very memory-efficient

### Pipeline takes too long

Estimated runtimes for a dataset of 10,000 sequences of ~1000bp:

| Method | Approximate time |
|---|---|
| `-d cosine` (default) | Minutes |
| `-d mash` | Minutes |
| `-d mmseqs2` | Minutes |
| `-d mmseqs-fast` | Seconds |
| `-d identity` | Hours |

If your run is too slow:
1. Switch to a faster distance method (`mmseqs-fast` is the fastest)
2. Increase threading (`-n 8` or higher)
3. Pre-reduce with Hobohm1 (`-a hobohm_split`)

## FAQ

**Q: What FASTA formats are supported?**
Standard multi-FASTA files with `>` headers. Both single-line and multi-line sequences are supported. The file can have any extension (`.fsa`, `.fasta`, `.fa`, `.fna`, `.faa`, etc.).

**Q: Can I use protein sequences?**
Yes, use `-s aminoacids`. Note that KMA-based distance methods (`cosine`, `jaccard`, etc.) only work with DNA. For proteins, use `-d mash`, `-d mmseqs2`, `-d mmseqs-fast`, or `-d identity`.

**Q: What is a good value for `-c`?**
It depends on your data. `-c 0.3` (70% identity) is a common starting point. If you're working with protein superfamilies, you might need `-c 0.5` or higher. If you're working with closely related gene variants, `-c 0.1` (90% identity) may be more appropriate. Run SpanSeq and check the bin balance in the stats file.

**Q: Can I reproduce a run?**
Yes. Given the same input FASTA, the same parameters, and the same tool versions, SpanSeq produces identical output. There is no randomness in the pipeline.

**Q: How many sequences can SpanSeq handle?**
It depends on the distance method:
- KMA, Mash, MMseqs2: comfortably handle tens of thousands of sequences
- MMseqs2 fast: handles hundreds of thousands
- GGSearch36: practical up to ~10,000 sequences (quadratic runtime)

**Q: What if I only need 2 bins (train/test)?**
Use `-b 2`. Or `-b 8,2` for an 80/20 split.

**Q: Can I add more sequences later?**
SpanSeq does not support incremental updates. You need to re-run the full pipeline with the expanded dataset.

## References

If using SpanSeq, please cite:

> Ferrer Florensa, Alfred, et al. "SpanSeq: Similarity-based sequence data splitting method for improved development and assessment of deep learning projects." *NAR Genomics and Bioinformatics*, Volume 6, Issue 3, September 2024. https://doi.org/10.1093/nargab/lqae106

SpanSeq relies on these external tools:

1. Clausen, Philip TLC, Frank M. Aarestrup, and Ole Lund. "Rapid and precise alignment of raw reads against redundant databases with KMA." *BMC bioinformatics* 19 (2018): 1-8.
2. Ondov, Brian D., et al. "Mash: fast genome and metagenome distance estimation using MinHash." *Genome biology* 17.1 (2016): 1-14.
3. Clausen, Philip TLC. "Scaling neighbor joining to one million taxa with dynamic and heuristic neighbor joining." *Bioinformatics* 39.1 (2023): btac774.
4. Hobohm, Uwe, et al. "Selection of representative protein data sets." *Protein Science* 1.3 (1992): 409-417.
5. Li, Weizhong, and Adam Godzik. "Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences." *Bioinformatics* 22.13 (2006): 1658-1659.
6. Steinegger, Martin, and Johannes Söding. "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets." *Nature biotechnology* 35.11 (2017): 1026-1028.
7. Pearson, William R. "Searching protein sequence libraries: comparison of the sensitivity and selectivity of the Smith-Waterman and FASTA algorithms." *Genomics* 11.3 (1991): 635-650.
