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
        nargs = "?", help = "sequence file either fasta or fastq")
parser.add_argument('-w', '--wanted_list', type = str,
        nargs = "?", help = "text file containing headers to be extracted")
parser.add_argument('-o', '--output', type = str,
        nargs = "?", help = "output file name")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.sequences is None or \
    args.wanted_list is None or \
    args.output is None:
    print("\n** required input missing\n"
          "** a sequence file, wanted list, and output name are all required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)


sequence_file = args.sequences
wanted_file = args.wanted_list
output_file = args.output

# get list wanted header names

wanted_list = [line.strip().lstrip(">").lstrip("@") for line in open(wanted_file) if line != ""]

# function to extract sequences in the wanted list for fasta or fastq

def extract_fasta_seqs(seq_fl, wanted, out_fl):
    count = 0
    fasta_seqs = SeqIO.parse(open(seq_fl), 'fasta')
    with open(out_fl, "w") as f:
        for seq in fasta_seqs:
            if seq.id in wanted:
                SeqIO.write([seq], f, "fasta")
                count += 1
    return(count)

def extract_fastq_seqs(seq_fl, wanted, out_fl):
    count = 0
    fastq_seqs = open(seq_fl)
    with open(out_fl, "w") as f:
        for title, seq, qual in FastqGeneralIterator(fastq_seqs):
            name = title.split()[0]
            if name in wanted:
                f.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                count += 1
    return(count)

# file type simply detected by looking at file ending

if sequence_file.endswith("fasta") or sequence_file.endswith("fa"):
    count = extract_fasta_seqs(sequence_file, wanted_list, output_file)
elif sequence_file.endswith("fastq") or sequence_file.endswith("fq"):
    count = extract_fastq_seqs(sequence_file, wanted_list, output_file)
else:
    print("** Error: Your sequence files must end with 'fasta' or 'fa' for fasta files, or 'fastq' or 'fq' for fastq "
          "files")
    sys.exit(1)

# print some output

print("\n** Saved {} records from {} to {}".format(count, sequence_file, output_file))
if count < len(wanted_list):
    print("** Warning: {} IDs not found in {}".format(len(wanted_list)-count, sequence_file))



