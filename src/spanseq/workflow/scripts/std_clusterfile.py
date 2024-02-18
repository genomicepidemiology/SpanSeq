import os
import pandas as pd
import csv
import subprocess
import argparse

def readlines_reversed(f):
    """ Iterate over the lines in a file in reverse. The file must be
    open in 'rb' mode. Yields the lines unencoded (as bytes), including the
    newline character. Produces the same result as readlines, but reversed.
    If this is used to reverse the line in a file twice, the result is
    exactly the same.
    """
    head = b""
    f.seek(0, 2)
    t = f.tell()
    buffersize, maxbuffersize = 64, 4096
    while True:
        if t <= 0:
            break
        # Read next block
        buffersize = min(buffersize * 2, maxbuffersize)
        tprev = t
        t = max(0, t - buffersize)
        f.seek(t)
        lines = f.read(tprev - t).splitlines(True)
        # Align to line breaks
        if not lines[-1].endswith((b"\n", b"\r")):
            lines[-1] += head  # current tail is previous head
        elif head == b"\n" and lines[-1].endswith(b"\r"):
            lines[-1] += head  # Keep \r\n together
        elif head:
            lines.append(head)
        head = lines.pop(0)  # can be '\n' (ok)
        # Iterate over current block in reverse
        for line in reversed(lines):
            yield line
    if head:
        yield head


def transform_clusterfiles(infile, outfile, fasta_file, dist):
    lines = read_file(infile)
    cluster = None
    n_neighbours = 0
    with open(outfile, "w") as output_f:
        descriptor = create_descr(fasta_file=fasta_file,
                            clstr_name=infile, dist=dist)
        output_f.write(descriptor)
        header = "#Sample\tNeighbors\tCluster\n"
        output_f.write(header)
        for l in lines:
            if l.startswith(">"):
                if cluster is not None:
                    clust_df = pd.DataFrame(sample_clst)
                    length_clus = len(clust_df)
                    clust_df[1] = length_clus-1
                    clust_df.to_csv(output_f, header=None, index=None,
                                    sep='\t', mode='a', quotechar=None,
                                    quoting=csv.QUOTE_NONE)
                cluster = l.replace(">Cluster ", "")
                sample_clst = []
            elif l[0].isdigit():
                end = l.split(", >")[1]
                fasta_header = end.split("... ")[0]
                sample_name = '"{}"'.format(fasta_header)
                sample_clst.append([sample_name, n_neighbours, cluster])
            else:
                raise TypeError("""CD-HIT file does not follow the correct"""
                                """ format""")
        clust_df = pd.DataFrame(sample_clst)
        length_clus = len(clust_df)
        clust_df[1] = length_clus-1
        clust_df.to_csv(output_f, header=None, index=None,
                        sep='\t', mode='a', quotechar=None,
                        quoting=csv.QUOTE_NONE)

def create_descr(fasta_file, clstr_name, dist):
    count_seq = 'grep ">" {} | wc -l'.format(fasta_file)
    num_seq = subprocess.check_output(count_seq, shell=True)
    with open(clstr_name, "rb") as clstr_file:
        f_rev = readlines_reversed(clstr_file)
        for line in f_rev:
            if line.startswith(b">"):
                num_cluster = line.rstrip(
                                    ).split(b">Cluster ")[-1].decode("utf-8")
                break
    dist_float = "{:.6f}".format(float(dist))

    string_descr = "# {0}\t{1}\t{2}\n".format(
                    num_seq.decode("utf-8").rstrip(), num_cluster, dist_float)
    return string_descr


def read_file(infile):
    with open(infile, "r") as input_f:
        line = True
        while line:
            line = input_f.readline().strip()
            if line:
                yield line

def arguments_std():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True,
                        help="Input file (output of cd-hit)")
    parser.add_argument("-o", "--output_file", required=True,
                        help="Output file (file format of ccphylo dbscan)")
    parser.add_argument("-f", "--fasta_file", required=True,
                        help="Fasta file (with sequences to cluster)")
    parser.add_argument("-d", "--dist_val", required=True,
                        help="Max distance value among clusters")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arguments_std()
    transform_clusterfiles(infile=args.input_file, outfile=args.output_file,
                           fasta_file=args.fasta_file, dist=args.dist_val)
