#!/usr/bin/env python
# python 3

# script to extract seqs from multi-fasta file
# and write to individual files
# Matthew J. Neave 2022-08-21

# library imports

import sys
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator # requires Biopython
from Bio import SeqIO

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("extract sequences from multi-fasta file into individual files\n"
                                 "note: ")

parser.add_argument('-s', '--sequences', type = str,
        nargs = "?", help = "sequence file either fasta or fastq")


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


# function to extract sequences in the wanted list for fasta or fastq

count = 0
fasta_seqs = SeqIO.parse(open(sequence_file), 'fasta')
for seq in fasta_seqs:
    seq.id = seq.id.replace("/", "_")
    with open(seq.id + ".fasta", "w") as f:
        SeqIO.write([seq], f, "fasta")
        count += 1


# print some output

print("\n** Saved {} records new files".format(count))


