import pandas as pd
import argparse
import sys
import os

def join_class(df_file, class_file, out_file):
    with open(df_file, "r") as dffile:
        header = dffile.readline().strip()
    df = pd.read_csv(df_file, sep="\t", skiprows=[0])
    class_df = pd.read_csv(class_file, sep="\t", header=None)
    change_cols = {}
    n_cols = len(list(class_df))
    for n in range(n_cols):
        old_n = list(class_df)[n]
        if n == 0:
            new_n = "Names"
        else:
            new_n = "Class_{}".format(n)
        change_cols[old_n] = new_n
    class_df.rename(columns=change_cols, inplace=True)
    final_df = pd.merge(df, class_df, how='inner', left_on="#Sample",
                        right_on="Names")
    del final_df["Names"]
    with open(out_file, 'w') as f:
        f.write(header)
        f.write("\n")
    if len(final_df) < len(df):
        raise ValueError("The file with classes contains less sequences than the original data")
    else:
        final_df.to_csv(out_file, index = False, sep="\t", mode='a')

def arguments_join():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_file", required=True,
                        help="Output folder (file format of ccphylo dbscan)")
    parser.add_argument("-d", "--df_file", required=True,
                        help="File with clusters")
    parser.add_argument("-c", "--classes_file", default=False,
                        help="File with classes")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arguments_join()
    join_class(df_file=args.df_file, class_file=args.classes_file,
                    out_file=args.output_file)

    