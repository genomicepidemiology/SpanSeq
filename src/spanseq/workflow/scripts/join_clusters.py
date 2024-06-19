import pandas as pd
import argparse
import sys
import os
from Bio import SeqIO


def create_fastas(dataframe, fasta_file, bins, outfolder):
    dataframe = pd.read_csv(dataframe, sep="\t")
    out_folder = os.path.abspath(outfolder)
    if os.path.isdir(out_folder):
        out_name = os.path.splitext(os.path.basename(fasta_file))[0]
        ext = "".join(os.path.splitext(os.path.basename(fasta_file))[1:])
    else:
        sys.exit("-o/--output_folder must be an existing folder")

    files_dict = {}
    for i in range(1,bins+1):
        out_bin = "{}_{}".format(out_name, i)
        out_path = "{path}/{name}_M{bin}{ext}".format(
                        path=out_folder, name=out_name,
                        bin=i, ext=ext)
        files_dict[out_bin] = {}
        files_dict[out_bin] = {"path":out_path, "exec":open(out_path, "w")}

    with open(fasta_file, "r") as fastafile:
        for record in SeqIO.parse(fastafile, "fasta"):
            partition = dataframe.loc[dataframe["id"]==record.id]["partition"].values
#            print(record.id, record.description, record)
            if len(partition) < 1:
                partition = dataframe.loc[dataframe["id"]==record.description]["partition"].values
            if len(partition) < 1:

                partition = dataframe[dataframe["id"].str.contains("{} ".format(record.description), regex=False)]["partition"].values
            if len(partition) > 1:
                partition = dataframe[dataframe["id"]==record.description]["partition"].values

            if len(partition) == 1:
                file_name = "{}_{}".format(out_name, int(partition))
                files_dict[file_name]["exec"].write(
                                            ">{}\n".format(record.description))
                files_dict[file_name]["exec"].write(
                                            "{}\n".format(record.seq))
            elif len(partition) > 1:
                sys.exit(("The fasta header %s is found more than once in the "
                         "cluster file" % record.description))
            else:
                print(record.id, record.description, record)
                continue
                sys.exit(("The fasta header %s is not found in the cluster file" % record.description))


def arguments_join():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_folder", required=True,
                        help="Output folder (file format of ccphylo dbscan)")
    parser.add_argument("-d", "--merged_df", required=True,
                        help="Merged dataframe dbscan+makespan")
    parser.add_argument("-f", "--fasta_file", default=False,
                        help="Fasta file (with sequences to cluster)")
    parser.add_argument("-b", "--bins", required=True, type=int,
                        help="Number of bins it has been divided the data in")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arguments_join()

    create_fastas(dataframe=args.merged_df, fasta_file=args.fasta_file, bins=args.bins,
                  outfolder=args.output_folder)
    
