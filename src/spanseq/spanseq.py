#    Copyright (c) 2023 Technical University of Denmark
#    All rights reserved.

 #   This source code is licensed under the BSD-style license found in the
 #   LICENSE file in the root directory of this source tree. 



import argparse
import subprocess
import os
import sys
from spanseq.utils.file_mixin import FileValidation
from spanseq.utils.environment_mixin import Environment
from spanseq.utils.configuration import Config
from spanseq.utils.pipeline_conf import Pipelines
from spanseq import __version__

# TODO: ADD clusterize          Cluster all thesequences, and distribute those clusters in the different bins/machines so there is an equal amount of clusters among machines.

class SpanSeq:

    AVAILABLE_ACTIONS = ["split", "reduce"]

    def __init__(self, safe_run=True):

        self.conda_env = Environment.get_envname()
        if safe_run:
            pass
            # TODO: Check if environment version is fine
            # Environment.check_version()

    @staticmethod
    def get_snakefilePath(file="workflow/Snakefile"):
        sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
        if not os.path.exists(sf):
            sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
        return sf

    @staticmethod
    def get_envfolderPath(folder="../../data/envs/"):
        sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), folder)
        if not os.path.isdir(sf):
            sys.exit("Unable to find the folder with environments; tried %s" % sf)
        return sf
    
    @staticmethod
    def get_configfilePath(file, folder=""):
        if folder == "":
            sf = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              file_loc)
        else:
            sf = folder
        if not os.path.isdir("{}/config".format(sf)):
            os.mkdir("{}/config".format(sf))
        file_loc = "{}/config/{}".format(sf, file)
        return file_loc


    @staticmethod
    def get_arguments():
        parser = argparse.ArgumentParser(prog="SpanSeq",
                                         usage='%(prog)s [options]',
                                         description='SpanSeq description', add_help=True)
        parent_parser = argparse.ArgumentParser(add_help=False)
        #General argumetns
        parser.add_argument("-v", "--version", action='version',
                            version=__version__)
        
        subparser = parser.add_subparsers(title="SPANSEQ MODES")
        ## General Data Arguments
        data_arg = parent_parser.add_argument_group('Common Data Options',
                            description="Data Options common with among all modes in SpanSeq. Include Input/Output options, as well as type of input or format of output")
        data_arg.add_argument("-i", "--input_fasta",
                            metavar="MULTIFASTA FILE",
                            type=lambda x: FileValidation.is_valid_file(file=x,
                                                                parser=splitter),
                            help=("Fasta file containing the dna or protein "
                            "sequences to divide/split in different machines."))
        data_arg.add_argument("-if", "--input_folder",
                            metavar="FOLDER FASTA",
                            type=FileValidation.readable_dir,
                            help=("Folder with fasta files containing the sequences to cluster."))
        data_arg.add_argument("-ib", "--input_batch",
                            metavar="LIST NAMES",
                            type=FileValidation.is_valid_file,
                            help=("File with the paths to the files containing the sequences to cluster"))
        data_arg.add_argument("-s", "--seqtype", type=str, required=True,
                            choices=["nucleotides", "aminoacids"], help=("Type of "
                            "sequence at the fasta file. The options are 'nucleotides'"
                            " and 'aminoacid'"))
        data_arg.add_argument("-o", "--output_folder", required=True,
                             help=("Output folder for temp files and result "
                             "files"))
        data_arg.add_argument("-f", "--output_format", default="minimal", type=str,
                             help=("Output format. 'minimal' returns the files prodcued by ccphylo to be mapped by the user (faster);"
                             " 'merged_table' returns a table with the name of the sequences and its partition and cluster assigned."
                             "'fasta_files' returns the table of 'merged_table' plus fasta files per partition (slower)."),
                             choices=["minimal", "merged_table", "fasta_files"])
        data_arg.add_argument("-tmp", "--temp_files",
                             help=("Location for temp files"), default=None)
        data_arg.add_argument("-r", "--keep_tmp", action="store_true",
                              help=("Do not remove temporary files"),
                              default=False)
        ## Software arguments
        softw_arg = parent_parser.add_argument_group('Common Software Options',
                            description="Options for overwriting software paths (instead of using Anaconda's installed)")
        softw_arg.add_argument("-KP", "--kmaPath", default=False, help=("Path to"
                              " the kma executable or the folder that contains"
                              " it. If used, overwrites the conda environment "
                              "that would be used instead. If used as a flag, "
                              "it will use the kma at the global path"))
        softw_arg.add_argument("-CP", "--ccphyloPath", default=False,
                              help=("Path to"
                              " the ccphylo executable or the folder that contains"
                              " it. If used as a flag, "
                              "it will use the ccphylo at the global path"))
        softw_arg.add_argument("-MP", "--mashPath", default=False, help=("Path to"
                              " the mash executable or the folder that contains"
                              " it. If used, overwrites the conda environment "
                              "that would be used instead. If used as a flag, "
                              "it will use the mash at the global path"))
        
        ## K-mer arguments
        kmer_arg = parent_parser.add_argument_group('Common K-mer Options',
                                    description="Common Options for K-mer indexing for KMA or Mash")
        kmer_arg.add_argument("-k", "--kmer_size", default=None, help=("Kmer "
                              "size used by MASH or KMA"), type=int)
        kmer_arg.add_argument("-m", "--minimizer_size", default=False, help=(
                              "Minimizer size used by kma"))   
        kmer_arg.add_argument("-p", "--prefix", default="-", help=("Use 'TG' when working with large databases (for KMA)"), type=str)
        kmer_arg.add_argument("-ME", "--MegaDB", action="store_true",
                              default=False,
                              help=("MEGADB. Use it for database of the order of 10⁶ (option for KMA)"))
        kmer_arg.add_argument("-l", "--max_length", type=int,
                              help=("Aproximate length of the longest sequence (option for Mash)"))

        ## Clustering arguments
        clustering_arg = parent_parser.add_argument_group('Common Distance and Partition arguments',
                                description="Minimum distance between partitions and amount of partitions created.")
        clustering_arg.add_argument("-c", "--min_dist", required=True, type=float,
                              help=("Maximum distance value between two"
                              " sequences allowed to be in different bins/machines. The "
                              "value must be between 0 and 1. The closer to one, the most restrictive the algorithm will be."))
        clustering_arg.add_argument("-b", "--bins", required=True,
                              help=("It can be a number or a list of numbers. "
                              "If it is a number, is the amount of "
                              "bins/machines the data is splitted on. If a "
                              "list, is the proportion of data in each "
                              "bin/machine"))
        ## Makespan arguments
        makespan_arg = parent_parser.add_argument_group('Makespan arguments')
        makespan_arg.add_argument("-mP", "--makespanProcess", choices=["DBF", "DFF"],
                              help=("Method on the makespan step. "
                              "The options can be DBF (Decreasing Best First/"
                              "Longest Processing Time (LPT)) or DFF ("
                              "Decreasing First Fit)"), default="DBF")
        makespan_arg.add_argument("-mW", "--makespanWeights",
                            choices=["none", "logX", "powX", "expX"],
                            help=("Weighting method on the makespan step. The "
                            "options can be:\nnone:	Do not weigh clusters\n"
                            "logX:	Weigh one plus logarithmicly with base X\n"
                            "powX:	Weigh polynomial with exponent X\n"
                            "expX:	Weigh exponential with exponential base X"),
                            default="none")
        ## Technical parameters for memory and speed consumption
        technical_arg = parent_parser.add_argument_group('Technical arguments')
        technical_arg.add_argument("-n", "--threads", default=1, help=("Threads used"
                              " by the pipeline"))
        ### SPLIT
        splitter = subparser.add_parser('split', help=("Distribute the "
                    "sequences in different bins/machines, so two sequences in different "
                    "bins are not similar above certain threshold. The amount of sequences"
                    " among the bins is as close as possible."), parents=[parent_parser])
        splitter.set_defaults(action="split")
        reduction = subparser.add_parser('reduce', help=("Remove sequences "
                        "from a database so no sequence is similar to any other sequence above"
                        " a certain threshold"), parents=[parent_parser])
        reduction.set_defaults(action="reduce")

        ## Pipeline arguments
        dist_arg = splitter.add_argument_group('Distance Arguments')
        dist_arg.add_argument("-a", "--approach", type=str,
                            choices=["all", "hobohm_reduce", "hobohm_split"],
                            default="all", help=("Mode in which SpanSeq split"
                            "is run. 'all' runs comparisions all vs.all, "
                            "while hobhom will run the hobhom algorithm."))
        dist_arg.add_argument("-d", "--distanceMethod",
                              default="szymkiewicz_simpson", type=str,
                              choices=["jaccard", "szymkiewicz_simpson",
                                    "cosine", "kmer_inv", "mash"],
                              help=("Method for calculating the distance "
                                    "between sequences by ccphylo"))
        dist_arg.add_argument("-hd", "--hobohm1_distance", default=False, help=( 
                            "Distance used for hobom1 algorithm to reduce the amount of sequence to be used for distance calculation."
                            "Only for --mode 'hobohm_reduce' and 'hobohm_split'. For example, if the value is 0.6, sequences with 0.6 "
                            "identical to another sequences will not be considered for distance calculation."), type=float)
        dist_arg.add_argument("-H", "--memory_disk", action="store_true",
                              help=("Allocate distance matrix on the disk"))


        ## Reduce options
        args = parser.parse_args()
        return args

    def set_configuration_run(self, config_filePath, args):
        configuration = Config()
        configuration.set_params(arg_params=args)

        configuration.set_environment()
        configuration.set_executables(args=args)

        configuration.write_config_file(out_file=config_filePath)
        return configuration

    def run_snakemake(self, args):
        if not os.path.isdir(args.output_folder):
            os.mkdir(args.output_folder)
        config_filePath = SpanSeq.get_configfilePath(file="config_spanseq.yaml",
                                             folder=args.output_folder)
        snakefilePath = SpanSeq.get_snakefilePath()
        envs_folder = SpanSeq.get_envfolderPath()

        configuration = self.set_configuration_run(
                                            config_filePath=config_filePath,
                                            args=args)

        log_folder = configuration.configfile["general"]["logdir"]
        cmd = ("snakemake --use-conda -p -c1 --snakefile {snakefile} "
               "--configfile {configfile} --conda-frontend mamba "
               "--conda-prefix {env_folder} --stats "
               "{statsfolder}snakemake.stats --directory {wrking_dir} "
               "--cores {threads}"
               ).format(snakefile=snakefilePath, configfile=config_filePath,
                    env_folder=envs_folder, statsfolder=log_folder,
                    wrking_dir=args.output_folder, threads=args.threads)
        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            print(e)
            exit(1)

def main():
    args = SpanSeq.get_arguments()
    run_spanseq = SpanSeq()
    if args.action in SpanSeq.AVAILABLE_ACTIONS:
        run_spanseq.run_snakemake(args=args)
    else:
        raise ValueError("The action {} is not available in SpanSeq at the moment".format(args.action))

if __name__ == '__main__':
    main()
