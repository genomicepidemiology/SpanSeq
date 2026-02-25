# SpanSeq

Homology-aware partitioning of biological sequence databases for unbiased machine learning.

SpanSeq splits a collection of DNA or protein sequences into bins such that the maximum similarity between any two sequences in different bins stays below a user-defined threshold. This prevents data leakage when building training, validation, and test sets.

## Installation

```bash
git clone https://github.com/genomicepidemiology/SpanSeq.git
cd SpanSeq

conda env create -n spanseq --file data/envs/spanseqenv.yml
conda activate spanseq

# cgecore v3 is not yet on PyPI — install directly from the development branch
pip install git+https://bitbucket.org/genomicepidemiology/cgecore.git@v3_dev

pip install -e .
```

## Quick Start

```bash
# Split nucleotide sequences into 3 bins with max 70% cross-bin similarity
spanseq split -i sequences.fsa -s nucleotides -o results/ -c 0.7 -b 3

# Get per-partition FASTA files
spanseq split -i sequences.fsa -s nucleotides -o results/ -c 0.7 -b 3 -f fasta_files

# Reduce redundancy
spanseq reduce -i sequences.fsa -s nucleotides -o results/ -c 0.7 -b 3
```

## Distance Methods

| Method | Tool | Sequence types | Flag |
|---|---|---|---|
| Cosine (default) | KMA | DNA | `-d cosine` |
| Jaccard | KMA | DNA | `-d jaccard` |
| Mash | Mash | DNA, protein | `-d mash` |
| Global identity | GGSearch36 | DNA, protein | `-d identity` |
| MMseqs2 identity | MMseqs2 | DNA, protein | `-d mmseqs2` |
| MMseqs2 fast cluster | MMseqs2 | DNA, protein | `-d mmseqs-fast` |

## Documentation

- [User Guide](docs/user_guide.md) — installation, CLI reference, examples
- [Developer Guide](docs/developer_guide.md) — architecture, adding distance methods, testing

## Citation

If using SpanSeq, please cite:

> Ferrer Florensa, Alfred, et al. "SpanSeq: Similarity-based sequence data splitting method for improved development and assessment of deep learning projects." *NAR Genomics and Bioinformatics*, Volume 6, Issue 3, September 2024. https://doi.org/10.1093/nargab/lqae106

## External Tools

SpanSeq uses:

1. [KMA](https://bitbucket.org/genomicepidemiology/kma) — Clausen et al. (2018)
2. [Mash](https://github.com/marbl/Mash) — Ondov et al. (2016)
3. [CCPhylo](https://bitbucket.org/genomicepidemiology/ccphylo) — Clausen (2023)
4. [CD-HIT](https://github.com/weizhongli/cdhit) — Li & Godzik (2006)
5. [GGSearch36](https://fasta.bioch.virginia.edu/) — Pearson (1991)
6. [MMseqs2](https://github.com/soedinglab/MMseqs2) — Steinegger & Söding (2017)

## License

Apache-2.0
