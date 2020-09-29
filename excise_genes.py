#!/usr/bin/env python
# python 3

# script to excise gene sequence from fasta file
# given gff list containing contig and start stop locations
# Matthew J. Neave 19.08.2020 <matthewjneave1@gmail.com>

# library imports

import sys
import argparse
from Bio import SeqIO # requires biopython

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("extract gene sequences from fasta file given gff list of wanted genes\n")

parser.add_argument('-s', '--sequences', type = str,
        nargs = "?", help = "sequence fasta file")
parser.add_argument('-g', '--gff_list', type = str,
        nargs = "?", help = "gff file containing genes to be extracted")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.sequences is None or \
    args.gff_list is None:
    print("\n** required input missing\n"
          "** a sequence file and gff list are all required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)


sequence_file = args.sequences
gff_file = args.gff_list

# get dict of wanted sequences

gff_dict = {}

with open(gff_file) as fl:
    for line in fl:
        if line.startswith("#"): continue
        line = line.strip()
        cols = line.split("\t")
        contig = cols[0]
        start = int(cols[3])
        stop = int(cols[4])
        product = cols[8].split(";")[-1].replace("product=", "")

        if contig in gff_dict:
            gff_dict[contig][start] = (start, stop, product)
        else:
            gff_dict[contig] = {start: (start, stop, product)}


# extract the sequences from fasta

count = 0

with open(sequence_file) as fl:
    for seq in SeqIO.parse(fl, "fasta"):
        if seq.id in gff_dict:
            for gene in gff_dict[seq.id]:

                start = gff_dict[seq.id][gene][0]
                end = gff_dict[seq.id][gene][1]
                product = gff_dict[seq.id][gene][2].replace(" ", "_").replace("/", "_").replace("'", "")

                fl_name = args.gff_list.split(".")[0]
                
                gene_seq = seq[start+1:end]
                gene_seq.id = fl_name + "_" + seq.id + "_" + str(start) + "_" + str(end) + "_" + product
                gene_seq.description = ""
                gene_output = open(fl_name + "_" + str(start) + "_" + product + ".fasta", "w")                
            
                SeqIO.write(gene_seq, gene_output, "fasta")
                count += 1


# print some output

print("\n** Saved {} records from {} to individual fasta files".format(count, sequence_file))


