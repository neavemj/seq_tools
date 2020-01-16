#!/usr/bin/env python
# python 3

# script to histogram of sequence lengths 
# Matthew J. Neave 17.01.2020

# library imports

import sys
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator # requires Biopython
from Bio import SeqIO
import numpy as np

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("draw histogram of sequence lengths from fasta or fastq file\n")

parser.add_argument('-s', '--sequences', type = str,
        help = "sequence file either fasta or fastq")
parser.add_argument('-b', '--bins', type = int, default=10,
        help = "sequence file either fasta or fastq")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.sequences is None:
    print("\n** required input missing\n"
          "** a sequence file is required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)


sequence_file = args.sequences

# function to get lengths for either fasta or fastq

def extract_fasta_seqs(seq_fl):
    lengths = []
    for record in SeqIO.parse(open(seq_fl), 'fasta'):
        lengths.append(len(record.seq))    
    return(lengths)

def extract_fastq_seqs(seq_fl):
    lengths = []
    for record in SeqIO.parse(open(seq_fl), 'fastq'):
        lengths.append(len(record.seq))    
    return(lengths)

# file type simply detected by looking at file ending

if sequence_file.endswith("fasta") or sequence_file.endswith("fa") or sequence_file.endswith("faa") or sequence_file.endswith("fna"):
    length_list = extract_fasta_seqs(sequence_file)
elif sequence_file.endswith("fastq") or sequence_file.endswith("fq"):
    length_list = extract_fastq_seqs(sequence_file)
else:
    print("** Error: Your sequence files must end with 'fasta' or 'fa' for fasta files, or 'fastq' or 'fq' for fastq "
          "files")
    sys.exit(1)

# use numpy to print histogram output
# this produces a couple of lists like this:
# (array([ 158,  266,  490, 1099, 2309, 3799, 4566, 4094, 2405, 1358]), array([1300., 1310., 1320., 1330., 1340., 1350., 1360., 1370., 1380., 1390., 1400.]))
# first list is the value, second is the outer bins (i.e. there is more bins than values)

hst = np.histogram(length_list, bins=args.bins)

for i, val in enumerate(hst[0]):
    bin1 = int(hst[1][i])
    bin2 = int(hst[1][i+1])
    # convert num of seqs into percentage of all seqs
    val_perc = round((val / sum(hst[0])) * 100) 
    val_stars = "*" * int(val_perc)
    # now nicely print this info off
    print("Bin {} to {}, {} ({}%):\n{}".format(bin1, bin2, val, val_perc, val_stars))

