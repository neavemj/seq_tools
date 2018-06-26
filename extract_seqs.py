#!/usr/bin/env python
# python 3

# script to extract seqs from fasta or fastq file
# given list of wanted headers
# Matthew J. Neave 26.06.2018 <matthewjneave1@gmail.com>

# library imports

import sys
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator # requires Biopython
from Bio import SeqIO

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("extract sequences from fasta or fastq file given list of sequence headers\n"
                                 "note: ")

parser.add_argument('-s', '--sequences', type = str,
        nargs = 1, help = "sequence file either fasta or fastq")
parser.add_argument('-w', '--wanted_list', type = str,
        nargs = 1, help = "text file containing headers to be extracted")
parser.add_argument('-o', '--output', type = str,
        nargs = 1, help = "output file name")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# get list and set of wanted header names

wanted_list = [line.strip().lstrip(">") for line in open(args.wanted_list[0]) if line != ""]
wanted_set = set(wanted_list)

# function to extract sequences in the wanted list

def extract_seqs(seq_fl, wanted, out_fl):
    count = 0
    fasta_seqs = SeqIO.parse(open(seq_fl), 'fasta')
    with open(out_fl, "w") as f:
        for seq in fasta_seqs:
            if seq.id in wanted:
                SeqIO.write([seq], f, "fasta")
                count += 1
    return(count)

count = extract_seqs(args.sequences[0], wanted_list, args.output[0])

print("Saved {} records from {} to {}".format(count, args.sequences[0], args.output[0]))
if count < len(wanted_list):
    print("Warning: {} IDs not found in {}".format(len(wanted_list)-count, args.sequences[0]))



