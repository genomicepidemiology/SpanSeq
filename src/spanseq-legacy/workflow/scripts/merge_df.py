import pandas as pd
import argparse


def join_clstr_bin(clstr_file, makespan_file):
    clstr_df = pd.read_csv(clstr_file, sep="\t", usecols=[0,2],
                           names=["id", "cluster"], skiprows=2,
                           dtype={"id": "string", "cluster": "int"})
    makespan_df = pd.read_csv(makespan_file, sep="\t", usecols=[0,3],
                              names=["cluster", "partition"], skiprows=1,
                              dtype={"cluster": "int", "partition": "int"})
    s1 = clstr_df.merge(makespan_df, left_on="cluster", right_on="cluster")

    return s1

def arguments_merge():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--clstr_file", required=True,
                        help="Input file (output of cd-hit)")
    parser.add_argument("-o", "--output_file", required=True,
                        help="Output folder (file format of ccphylo dbscan)")
    parser.add_argument("-m", "--makespan_file", required=True,
                        help="File produced by makespan")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arguments_merge()
    df_join = join_clstr_bin(clstr_file=args.clstr_file,
                             makespan_file=args.makespan_file)
    df_join.to_csv(args.output_file, sep="\t", index=False)