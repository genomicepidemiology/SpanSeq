SpanSeq
==============

Spanseq is an easy-to-use software to reduce the effect of homology in databases of biological sequences for the development of machine learning methods.


# Setup

## Install Anaconda

In order to set up SpanSeq, Anaconda is required. All the other dependencies can be installed manually or through Anaconda. SpanSeq is based on Snakemake which allows to run steps of the workflow in parallel on a cluster.

You need to install anaconda or miniconda. If you haven't done it already you need to configure conda with the bioconda-channel and the conda-forge channel. This are soruces for packages beyond the default one.

```console
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

## Install SpanSeq

### Install SpanSeq from Bitbucket

To run SpanSeq locally, this repository has to be cloned, a bioconda environment has to be activated and the package installed.

```console
git clone https://bitbucket.org/genomicepidemiology/spanseq.git
cd spanseq

conda env create -n spanseq --file data/envs/spanseqenv.yml

pip install -e .
```
### Install SpanSeq from Anaconda

## Use test

To test the installation, run

```console
conda activate spanseq

spanseq split -d 0.7 -i ./test/data/aminoglycoside.fsa -p 1 -s dna -o ./data/ -b 3
```
# Getting Started

SpanSeq work with two modes: split and reduction

## SpanSeq Split

SpanSeq Split divides the database of dna or protein sequences in a number of machines or bins, so the maximum homology between two sequences in different bins is never above a certain threshold.

### Command Options
```console
usage: SpanSeq [options] split [-h] [-i MULTIFASTA FILE] [-if FOLDER FASTA] [-ib LIST NAMES] -s {nucleotides,aminoacids} -o OUTPUT_FOLDER [-f {minimal,merged_table,fasta_files}] [-tmp TEMP_FILES] [-r]
                               [-KP KMAPATH] [-CP CCPHYLOPATH] [-MP MASHPATH] [-k KMER_SIZE] [-m MINIMIZER_SIZE] [-p PREFIX] [-ME] [-l MAX_LENGTH] -c MIN_DIST -b BINS [-mP {DBF,DFF}]
                               [-mW {none,logX,powX,expX}] [-n THREADS] [-a {all,hobohm_reduce,hobohm_split}] [-d {jaccard,szymkiewicz_simpson,cosine,kmer_inv,mash}] [-hd HOBOHM1_DISTANCE] [-H]

options:
  -h, --help            show this help message and exit

Common Data Options:
  Data Options common with among all modes in SpanSeq. Include Input/Output options, as well as type of input or format of output

  -i MULTIFASTA FILE, --input_fasta MULTIFASTA FILE
                        Fasta file containing the dna or protein sequences to divide/split in different machines.
  -if FOLDER FASTA, --input_folder FOLDER FASTA
                        Folder with fasta files containing the sequences to cluster.
  -ib LIST NAMES, --input_batch LIST NAMES
                        File with the paths to the files containing the sequences to cluster
  -s {nucleotides,aminoacids}, --seqtype {nucleotides,aminoacids}
                        Type of sequence at the fasta file. The options are 'nucleotides' and 'aminoacid'
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder for temp files and result files
  -f {minimal,merged_table,fasta_files}, --output_format {minimal,merged_table,fasta_files}
                        Output format. 'minimal' returns the files prodcued by ccphylo to be mapped by the user (faster); 'merged_table' returns a table with the name of the sequences and its partition
                        and cluster assigned.'fasta_files' returns the table of 'merged_table' plus fasta files per partition (slower).
  -tmp TEMP_FILES, --temp_files TEMP_FILES
                        Location for temp files
  -r, --keep_tmp        Do not remove temporary files

Common Software Options:
  Options for overwriting software paths (instead of using Anaconda's installed)

  -KP KMAPATH, --kmaPath KMAPATH
                        Path to the kma executable or the folder that contains it. If used, overwrites the conda environment that would be used instead. If used as a flag, it will use the kma at the
                        global path
  -CP CCPHYLOPATH, --ccphyloPath CCPHYLOPATH
                        Path to the ccphylo executable or the folder that contains it. If used as a flag, it will use the ccphylo at the global path
  -MP MASHPATH, --mashPath MASHPATH
                        Path to the mash executable or the folder that contains it. If used, overwrites the conda environment that would be used instead. If used as a flag, it will use the mash at the
                        global path

Common K-mer Options:
  Common Options for K-mer indexing for KMA or Mash

  -k KMER_SIZE, --kmer_size KMER_SIZE
                        Kmer size used by MASH or KMA
  -m MINIMIZER_SIZE, --minimizer_size MINIMIZER_SIZE
                        Minimizer size used by kma
  -p PREFIX, --prefix PREFIX
                        Use 'TG' when working with large databases (for KMA)
  -ME, --MegaDB         MEGADB. Use it for database of the order of 10⁶ (option for KMA)
  -l MAX_LENGTH, --max_length MAX_LENGTH
                        Aproximate length of the longest sequence (option for Mash)

Common Distance and Partition arguments:
  Minimum distance between partitions and amount of partitions created.

  -c MIN_DIST, --min_dist MIN_DIST
                        Maximum distance value between two sequences allowed to be in different bins/machines. The value must be between 0 and 1. The closer to one, the most restrictive the algorithm will be.
  -b BINS, --bins BINS  It can be a number or a list of numbers. If it is a number, is the amount of bins/machines the data is splitted on. If a list, is the proportion of data in each bin/machine

Makespan arguments:
  -mP {DBF,DFF}, --makespanProcess {DBF,DFF}
                        Method on the makespan step. The options can be DBF (Decreasing Best First/Longest Processing Time (LPT)) or DFF (Decreasing First Fit)
  -mW {none,logX,powX,expX}, --makespanWeights {none,logX,powX,expX}
                        Weighting method on the makespan step. The options can be: none: Do not weigh clusters logX: Weigh one plus logarithmicly with base X powX: Weigh polynomial with exponent X expX:
                        Weigh exponential with exponential base X

Technical arguments:
  -n THREADS, --threads THREADS
                        Threads used by the pipeline

Distance Arguments:
  -a {all,hobohm_reduce,hobohm_split}, --approach {all,hobohm_reduce,hobohm_split}
                        Mode in which SpanSeq splitis run. 'all' runs comparisions all vs.all, while hobhom will run the hobhom algorithm.
  -d {jaccard,szymkiewicz_simpson,cosine,kmer_inv,mash}, --distanceMethod {jaccard,szymkiewicz_simpson,cosine,kmer_inv,mash}
                        Method for calculating the distance between sequences by ccphylo
  -hd HOBOHM1_DISTANCE, --hobohm1_distance HOBOHM1_DISTANCE
                        Distance used for hobom1 algorithm to reduce the amount of sequence to be used for distance calculation.Only for --mode 'hobohm_reduce' and 'hobohm_split'. For example, if the value is 0.6, sequences with 0.6 identical to another sequences will not be considered for distance calculation. Notice that as higher this value, the most restrictive the clustering will be.
  -H, --memory_disk     Allocate distance matrix on the disk

```
### Approaches
SpanSeq has three distance calculation approaches:
#### ALL
All the input sequences will be considered when building the distance matrix.
#### Hobohm_reduce, Hobohm_split
Previous to creating the distance matrix, a Hobohm 1 clustering[5] step is applied to the sequences, and the representatives of those clusters are used for distance calculation. When using that option, it is necessary to introduce a value for -hd/--hobohm1_distance, which should be higher than the distance for clustering later. If using hobohm_reduce, the only the representatives of the Hobohm 1 Clustering step will be reported, but if hobohm_split, the clustered sequences around representatives will be later merged with the clusters created with DBScan. This step is only available on nucleotide sequences.

### Pipelines

SpanSeq can use a personally made pipeline, but it is recommended to use one of the pre-made pipelines. Below there is a detailed descirption of them.

#### Scheme 1

This pipeline can work for DNA sequences. It works running the methods below in this order:
  1.	`kma index` - Index the Fasta File with KMA[2]. If Hobohm 1 clustering is used, the step is performed while indexing.
  2.	`kma dist` - Calculates distances between templates from kma index. The standard distance method is Jaccard distance.
  3.	`ccphylo dbscan` - Make a DBSCAN given a set of phylip distance matrices using CCPhylo[3], returning the data in different clusters. The parameter of 'minimum neighbour' is set to 1.
  4.	`ccphylo makespan` - Distribute clusters with different amounts of sequences in a certain amount of bins or machines, so the amount of sequences among bins/machines is as close as possible to equal.

#### Scheme 2
This pipeline can work for DNA and protein sequences. It works running the methods below in this order:
  1.	`mash` - Calculates distances between templates from fasta file using Mash[4]. The standard distance method is Jaccard distance.
  2.	`ccphylo dbscan` - Make a DBSCAN given a phylip distance matrix, returning the data in different clusters. The parameter of 'minimum neighbour' is set to 1.
  3.	`ccphylo makespan` - Distribute clusters with different amounts of sequences in a certain amount of bins or machines, so the amount of sequences among bins/machines is as close as possible to equal.


## SpanSeq Reduce (Beta)
SpanSeq reduces the database of dna or protein sequences in a number of machines or bins, so the maximum homology between two sequences in different bins is never above a certain threshold.

### Command Options
```console
usage: SpanSeq [options] reduce [-h] [-i MULTIFASTA FILE] [-if FOLDER FASTA] [-ib LIST NAMES] -s {nucleotides,aminoacids} -o OUTPUT_FOLDER [-f {minimal,merged_table,fasta_files}] [-tmp TEMP_FILES] [-r]
                                [-KP KMAPATH] [-CP CCPHYLOPATH] [-MP MASHPATH] [-k KMER_SIZE] [-m MINIMIZER_SIZE] [-p PREFIX] [-ME] [-l MAX_LENGTH] -c MIN_DIST -b BINS [-mP {DBF,DFF}]
                                [-mW {none,logX,powX,expX}] [-n THREADS]

options:
  -h, --help            show this help message and exit

Common Data Options:
  Data Options common with among all modes in SpanSeq. Include Input/Output options, as well as type of input or format of output

  -i MULTIFASTA FILE, --input_fasta MULTIFASTA FILE
                        Fasta file containing the dna or protein sequences to divide/split in different machines.
  -if FOLDER FASTA, --input_folder FOLDER FASTA
                        Folder with fasta files containing the sequences to cluster.
  -ib LIST NAMES, --input_batch LIST NAMES
                        File with the paths to the files containing the sequences to cluster
  -s {nucleotides,aminoacids}, --seqtype {nucleotides,aminoacids}
                        Type of sequence at the fasta file. The options are 'nucleotides' and 'aminoacid'
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder for temp files and result files
  -f {minimal,merged_table,fasta_files}, --output_format {minimal,merged_table,fasta_files}
                        Output format. 'minimal' returns the files prodcued by ccphylo to be mapped by the user (faster); 'merged_table' returns a table with the name of the sequences and its partition
                        and cluster assigned.'fasta_files' returns the table of 'merged_table' plus fasta files per partition (slower).
  -tmp TEMP_FILES, --temp_files TEMP_FILES
                        Location for temp files
  -r, --keep_tmp        Do not remove temporary files

Common Software Options:
  Options for overwriting software paths (instead of using Anaconda's installed)

  -KP KMAPATH, --kmaPath KMAPATH
                        Path to the kma executable or the folder that contains it. If used, overwrites the conda environment that would be used instead. If used as a flag, it will use the kma at the
                        global path
  -CP CCPHYLOPATH, --ccphyloPath CCPHYLOPATH
                        Path to the ccphylo executable or the folder that contains it. If used as a flag, it will use the ccphylo at the global path
  -MP MASHPATH, --mashPath MASHPATH
                        Path to the mash executable or the folder that contains it. If used, overwrites the conda environment that would be used instead. If used as a flag, it will use the mash at the
                        global path

Common K-mer Options:
  Common Options for K-mer indexing for KMA or Mash

  -k KMER_SIZE, --kmer_size KMER_SIZE
                        Kmer size used by MASH or KMA
  -m MINIMIZER_SIZE, --minimizer_size MINIMIZER_SIZE
                        Minimizer size used by kma
  -p PREFIX, --prefix PREFIX
                        Use 'TG' when working with large databases (for KMA)
  -ME, --MegaDB         MEGADB. Use it for database of the order of 10⁶ (option for KMA)
  -l MAX_LENGTH, --max_length MAX_LENGTH
                        Aproximate length of the longest sequence (option for Mash)

Common Distance and Partition arguments:
  Minimum distance between partitions and amount of partitions created.

  -c MIN_DIST, --min_dist MIN_DIST
                        Maximum distance value between two sequences allowed to be in different bins/machines. The value must be between 0 and 1. The closer to one, the most restrictive the algorithm
                        will be.
  -b BINS, --bins BINS  It can be a number or a list of numbers. If it is a number, is the amount of bins/machines the data is splitted on. If a list, is the proportion of data in each bin/machine

Makespan arguments:
  -mP {DBF,DFF}, --makespanProcess {DBF,DFF}
                        Method on the makespan step. The options can be DBF (Decreasing Best First/Longest Processing Time (LPT)) or DFF (Decreasing First Fit)
  -mW {none,logX,powX,expX}, --makespanWeights {none,logX,powX,expX}
                        Weighting method on the makespan step. The options can be: none: Do not weigh clusters logX: Weigh one plus logarithmicly with base X powX: Weigh polynomial with exponent X expX:
                        Weigh exponential with exponential base X

Technical arguments:
  -n THREADS, --threads THREADS
                        Threads used by the pipeline
```

# References
If using SpanSeq, please cite:

Spanseq uses the other softwares and algorithms:
2. Clausen, Philip TLC, Frank M. Aarestrup, and Ole Lund. "Rapid and precise alignment of raw reads against redundant databases with KMA." BMC bioinformatics 19 (2018): 1-8.
3. Ondov, Brian D., et al. "Mash: fast genome and metagenome distance estimation using MinHash." Genome biology 17.1 (2016): 1-14.
4. Clausen, Philip TLC. "Scaling neighbor joining to one million taxa with dynamic and heuristic neighbor joining." Bioinformatics 39.1 (2023): btac774.
5. Hobohm, Uwe, et al. "Selection of representative protein data sets." Protein Science 1.3 (1992): 409-417.
