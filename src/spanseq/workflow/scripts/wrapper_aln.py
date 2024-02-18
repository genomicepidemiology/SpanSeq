import os
import numpy as np
import subprocess
import sys
import json
import random
import time
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class _Aligners:

    pass


class GGSearch36_Aligner(_Aligners):

    @staticmethod
    def run_ggsearch(input, target, outfile, max_len, cores, seq_type,
                     evalue=20000, loc_ggsearch=""):
        if seq_type == "nucleotides":
            seq_str = "-n"
        elif seq_type == "aminoacids":
            seq_str = "-p"
        else:
            raise ValueError("Seq_type can only be 'dna' or 'protein'")

        command = ('{loc_ggsearch}ggsearch36 "{input}" "{target}" {seq_str} '
                   '-m 9i -T {cores} -E {evalue} -M 1-{max_len} -d 0 > '
                   '"{outfile}"'.format(
                        loc_ggsearch=loc_ggsearch, input=input, target=target,
                        outfile=outfile, max_len=max_len, seq_str=seq_str,
                        cores=cores, evalue=evalue))
        run = subprocess.Popen(command, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, shell=True)
        stdout, stderr = run.communicate()
        p_status = run.wait()

    @staticmethod
    def read_ggsearchfile(ggsearch_file, score_bool=False, len_bool=False):
        aln_metadata = {}
        seq1 = str(os.path.basename(ggsearch_file).replace(".ggsearch36", ""))
        if os.path.isfile(ggsearch_file):
            with open(ggsearch_file, "r") as ggsearchfile:
                for line in ggsearchfile:
                    if line.startswith("The best scores are:"):
                        while True:
                            new_line = next(ggsearchfile)
                            if new_line.startswith(">>>"):
                                break
                            info_aln = new_line.strip().split(" ")
                            subject = info_aln[0]
                            id = ""
                            n_wscore = ""
                            pos = 0
                            for n in reversed(info_aln):
                                if n == "":
                                    continue
                                else:
                                    if pos == 2:
                                        id = float(n)
                                    elif pos == 5:
                                        n_wscore = float(n)
                                    elif pos == 0:
                                        len_aln = float(n)
                                    pos += 1

                            if subject in aln_metadata:
                                if id > aln_metadata[subject]["identity"]:
                                    aln_metadata[subject]["identity"] = id
                                    if score_bool:
                                        aln_metadata[subject]["score"] = n_wscore
                                    if len_bool:
                                        aln_metadata[subject]["length"] = len_aln

                                else:
                                    continue
                            else:
                                aln_metadata[subject] = {}
                                aln_metadata[subject]["identity"] = id
                                if score_bool:
                                    aln_metadata[subject]["score"] = n_wscore
                                if len_bool:
                                    aln_metadata[subject]["length"] = len_aln

        seq_metadata = {seq1: aln_metadata}
        return seq_metadata


class Aligner_Wrapper:

    # TODO: Seq name longer than 30 charac in emboss

    aligners_avail = ["ggsearch36", "emboss"]
    seqtypes_avail = ["nucleotides", "aminoacids"]

    def __init__(self, aligner, out_folder, seq_type, loc_aligner=""):

        if aligner in Aligner_Wrapper.aligners_avail:
            self.aligner = aligner
        else:
            raise ValueError("The aligner {} is not available in SpanSeq. The \
                              current available aligners are: {}".format(
                              aligner, ", ".join(
                              Aligner_Wrapper.aligners_avial)))
        if seq_type in Aligner_Wrapper.seqtypes_avail:
            self.seq_type = seq_type
        else:
            raise ValueError("The sequence type {} is not available in \
                              Spanseq. The current available sequence types \
                              are: 'nucleotides' and 'aminoacids'.")
        self.location_aln = loc_aligner
        self.results_folder = os.path.join(
                                    out_folder, "{}_results".format(aligner))
        if not os.path.isdir(self.results_folder):
            os.mkdir(self.results_folder)
        self.out_folder = out_folder

    @staticmethod
    def concatenate_files(target_file, input_file):
        command = 'cat "{}" >> {}'.format(input_file, target_file)
        run = subprocess.Popen(command, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, shell=True)
        stdout, stderr = run.communicate()
        p_status = run.wait()

    @staticmethod
    def create_targetfile(tmp_folder):
        target_file = os.path.relpath(os.path.join(tmp_folder,
                                                            "target_file.fsa"))
        open(target_file, 'w').close()
        return target_file

    @staticmethod
    def print_counter(i):
        sys.stdout.write("\rSequences processed %i" % i)
        sys.stdout.flush()

    @staticmethod
    def rewrite_firstline(file, number):
        print(number)
        new_line = "\t{}".format(number)
        cmd = ['sed', '-i', '-e', '1,1s/.*/' + new_line + '/g', str(file)]

        subprocess.call(cmd)

    @staticmethod
    def write_matrix(idenmatrix_file, aln_data, seq_names, inverse=False,
                     scorematrix_file=False, lenmatrix_file=False,
                     far_fill={"identity":-1, "score":"nan", "aln_len": "nan"}):
        for name in aln_data:
            seq1 = name.replace("/", "-").strip()
            idenmatrix_file.write("{}".format(seq1))
            if scorematrix_file:
                scorematrix_file.write("{}".format(seq1))
            if lenmatrix_file:
                lenmatrix_file.write("{}".format(seq1))
            if aln_data[name]:
                for seq2 in seq_names:
                    if len(seq2) < 30:
                        try:
                            datasub = aln_data[seq1][seq2]
                        except KeyError:
                            datasub = far_fill
                    else:
                        try:
                            datasub = aln_data[seq1][seq2[:31]]
                        except KeyError:
                            datasub = far_fill
                    if inverse:
                        identity_data = float(datasub["identity"])
                    else:
                        identity_data = 1-float(datasub["identity"])
                    idenmatrix_file.write("\t{}".format(identity_data))

                    if scorematrix_file:
                        score_data = float(datasub["score"])
                        scorematrix_file.write("\t{}".format(score_data))

                    if lenmatrix_file:
                        len_data = float(datasub["length"])
                        lenmatrix_file.write("\t{}".format(len_data))

            idenmatrix_file.write("\n")
            if scorematrix_file:
                scorematrix_file.write("\n")
            if lenmatrix_file:
                lenmatrix_file.write("\n")

    def run_aligner(self, fasta_file, phy_iden, phy_score=False,
                    phy_length=False, tmp_folder=None, test=False, **kwargs):
        if tmp_folder is None:
            tmp_folder = os.path.join(self.out_folder, "tmp_aligner")
        target_file = Aligner_Wrapper.create_targetfile(tmp_folder=tmp_folder)
        seen_sequences = []
        count = 0

        idenmatrix = open(phy_iden, "w")
        idenmatrix.write("\t{}\n".format("0"))

        if phy_score:
            scorematrix = open(phy_score, "w")
            scorematrix.write("\t{}\n".format("0"))
            phy_score = scorematrix
        if phy_length:
            lenmatrix = open(phy_length, "w")
            lenmatrix.write("\t{}\n".format("0"))
            phy_length = lenmatrix

        for record in SeqIO.parse(fasta_file, "fasta"):
            writteable_record = str(record.id).replace("/", "-")
            if test:
                fasta_tmpfile = os.path.relpath("{}/{}.fsa".format(tmp_folder,
                                                            writteable_record))
            else:
                fasta_tmpfile = os.path.relpath("{}/in_sq.fsa".format(
                                                                tmp_folder))

            with open(fasta_tmpfile, "w") as output_handle:
                SeqIO.write(record, output_handle, "fasta")
            result_file = os.path.relpath("{}/{}.{}".format(
                                self.results_folder, writteable_record,
                                self.aligner))

            if self.aligner == "emboss":
                Emboss_Aligner.run_emboss(input=fasta_tmpfile,
                                target=target_file, outfile=result_file,
                                seq_type=self.seq_type,
                                loc_emboss=self.location_aln)
                aln_data = Emboss_Aligner.read_embossfile(
                                                emboss_file=result_file)
            elif self.aligner == "ggsearch36":
                GGSearch36_Aligner.run_ggsearch(input=fasta_tmpfile,
                            target=target_file, outfile=result_file,
                            seq_type=self.seq_type,
                            loc_ggsearch=self.location_aln, **kwargs)
                aln_data = GGSearch36_Aligner.read_ggsearchfile(
                                            ggsearch_file=result_file,
                                            score_bool=phy_score,
                                            len_bool=phy_length)
            else:
                sys.exit("No aligner with that name")
            Aligner_Wrapper.concatenate_files(target_file=target_file,
                                              input_file=fasta_tmpfile)
            Aligner_Wrapper.write_matrix(idenmatrix_file=idenmatrix,
                                        scorematrix_file=phy_score,
                                        lenmatrix_file=phy_length,
                                        aln_data=aln_data,
                                        seq_names=seen_sequences)
            seen_sequences.append(writteable_record)
            if not test:
                os.remove(fasta_tmpfile)
                os.remove(result_file)
            count += 1
            Aligner_Wrapper.print_counter(count)



        idenmatrix.close()
        Aligner_Wrapper.rewrite_firstline(file=phy_iden, number=count)
        if phy_score:
            scorematrix.close()
            Aligner_Wrapper.rewrite_firstline(file=phy_score, number=count)

        if not test:
            os.remove(target_file)
            os.rmdir(tmp_folder)
        sys.stdout.write("\n")


    @staticmethod
    def set_arguments():
        parser = argparse.ArgumentParser()
        parser.add_argument("-a", "--aligner", required=True,
                            help="Aligner software used",
                            choices={"emboss", "ggsearch36"})
        parser.add_argument("-o", "--output_folder", required=True,
                            help="Output folder")
        parser.add_argument("-t", "--tmp_folder", help="Temporary folder")
        parser.add_argument("-s", "--seq_type", required=True,
                            help="Sequence type",
                            choices={"nucleotides", "aminoacids"})
        parser.add_argument("-f", "--fasta_file", required=True,
                            help="Fasta file (with sequences to align)")
        parser.add_argument("-S", "--score", action="store_true",
                            help=("Output also a score matrix, besides the \
                                   identity matrix"))
        parser.add_argument("-L", "--aln_len", action="store_true",
                            help=("Output also a alignment length matrix, \
                                    besides the identity matrix"))
        parser.add_argument("-I", "--inverse", action="store_true",
                            help=("When used, the matrix will be done with the"
                            " inverse (1-) of the identities"))
        parser.add_argument("-l", "--aligner_loc", help=("Folder where the "
                            "aligner is located"), type=str, default="")
        parser.add_argument("-m", "--max_len", help=("Maximum length of a "
                            "sequence to be aligned (only for ggsearch36)"),
                            type=int)
        parser.add_argument("-p", "--cores", type=int, default=None,
                            help=("Cores to be used (Only for ggsearch36). "
                                  "By default will be set to 1"))
        parser.add_argument("-e", "--evalue", type=float, default=None,
                            help=("Max e-value to report a hit "
                                  "(only for ggsearch36). By default will be "
                                  "set to 20000"))
        args = parser.parse_args()
        args_fixed = Aligner_Wrapper.config_arguments(args)
        return args_fixed

    @staticmethod
    def config_arguments(args):
        config_aligner = {}
        config_aligner["seq_type"] = args.seq_type
        if not os.path.isdir(args.output_folder):
            os.mkdir(args.output_folder)
        config_aligner["output_folder"] = args.output_folder


        if args.tmp_folder is None:
            args.tmp_folder = os.path.join(args.output_folder, "tmp")
        if not os.path.isdir(args.tmp_folder):
            os.mkdir(args.tmp_folder)
        config_aligner["tmp_folder"] = args.tmp_folder

        if not os.path.isfile(args.fasta_file):
            raise OSError("The file {} does not exists".format(args.fasta_file))

        config_aligner["input_file"] = args.fasta_file
        name_file = os.path.basename(args.fasta_file).split(".")[0]
        config_aligner["matrix_iden"] = os.path.join(args.output_folder,
                                                     name_file + ".phy")
        if args.score:
            config_aligner["matrix_score"] = os.path.join(args.output_folder,
                                                     name_file + "_score.phy")
        else:
            config_aligner["matrix_score"] = args.score

        if args.aln_len:
            config_aligner["matrix_alnlen"] = os.path.join(args.output_folder,
                                                     name_file + "_alnlen.phy")
        else:
            config_aligner["matrix_alnlen"] = args.aln_len

        config_aligner["loc_aligner"] = args.aligner_loc

        if args.aligner == "emboss":
            config_aligner["aligner"] = "emboss"
            config_aligner["params"] = {}

            if args.max_len is not None:
                raise ValueError("The option -m/--max_len is only available "
                                 "for ggsearch36")
            if args.cores is not None:
                raise ValueError("The option -p/--cores is only available "
                                 "for ggsearch36")
            if args.evalue is not None:
                raise ValueError("The option -e/--evalue is only available "
                                 "for ggsearch36")
        elif args.aligner == "ggsearch36":
            config_aligner["aligner"] = "ggsearch36"
            config_aligner["params"] = {}
            if args.max_len is None:
                raise ValueError("The option -m/--max_len needs to be "
                                 "specified when running ggsearch36")
            config_aligner["params"]["max_len"] = args.max_len
            if args.cores is None:
                args.cores = 1
            config_aligner["params"]["cores"] = args.cores

            if args.evalue is None:
                args.evalue = 20000
            config_aligner["params"]["evalue"] = args.evalue

        if args.inverse:
            config_aligner["params"]["inverse"] = True
        else:
            config_aligner["params"]["inverse"] = False
        return config_aligner


if __name__ == '__main__':
    config_aligner = Aligner_Wrapper.set_arguments()
    aligner_run = Aligner_Wrapper(aligner=config_aligner["aligner"],
                                  out_folder=config_aligner["output_folder"],
                                  seq_type=config_aligner["seq_type"],
                                  loc_aligner=config_aligner["loc_aligner"])
    if config_aligner["aligner"] == "emboss":
        aligner_run.run_aligner(fasta_file=config_aligner["input_file"],
                                phy_iden=config_aligner["matrix_iden"],
                                phy_score=config_aligner["matrix_score"],
                                tmp_folder=config_aligner["tmp_folder"])
    elif config_aligner["aligner"] == "ggsearch36":
        aligner_run.run_aligner(fasta_file=config_aligner["input_file"],
                                phy_iden=config_aligner["matrix_iden"],
                                phy_score=config_aligner["matrix_score"],
                                phy_length=config_aligner["matrix_alnlen"],
                                tmp_folder=config_aligner["tmp_folder"],
                                max_len=config_aligner["params"]["max_len"],
                                cores=config_aligner["params"]["cores"],
                                evalue=config_aligner["params"]["evalue"])
