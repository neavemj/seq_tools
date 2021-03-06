#!/usr/bin/env python
# python 3

# script to calculate and plot the gc distribution of a fasta or fastq file
# Matthew J. Neave 20.07.2018 <matthewjneave1@gmail.com>

# library imports

import sys
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
import pandas as pd
import matplotlib.pyplot as plt

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("calculate and plot gc distribution of fasta or fastq file\n")

parser.add_argument('-s', '--sequence_file', type = str,
                    nargs = "?", help = "fasta or fastq file")
parser.add_argument('-o', '--output_prefix', type = str,
                    nargs = "?", help = "a name to prefix the output files")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.sequence_file is None or args.output_prefix is None:
    print("\n** required input missing\n"
          "** a sam file and output prefix are all required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

seq_fl = args.sequence_file
output_prefix = args.output_prefix

# function for a fasta file

def fasta_gc(fasta):
    gc_dict = {}
    for record in SeqIO.parse(fasta, "fasta"):
        gc = round(GC(record.seq), 1)
        if gc in gc_dict:
            gc_dict[gc] += 1
        else:
            gc_dict[gc] = 1
    return(gc_dict)


gc_dict = fasta_gc(seq_fl)
gc_values = list(gc_dict.keys())
gc_count = list(gc_dict.values())

gc_df = pd.DataFrame({'gc': gc_values, 'gc_count': gc_count})
gc_df = gc_df.sort_values(by=['gc'])

print(gc_df)

gc_df.to_csv(output_prefix + ".gc_hist.txt", sep="\t", index=False)

gc_df.plot(x="gc", y="gc_count", kind="line", xlim=[0,100])
plt.show()
