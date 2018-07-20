#!/usr/bin/env python
# python 3

# script to convert a sam file to fastq
# Matthew J. Neave 20.07.2018 <matthewjneave1@gmail.com>

# library imports

import sys
import argparse

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("convert sam file to forward, reverse and single fastq files\n")

parser.add_argument('-s', '--sam_file', type = str,
                    nargs = "?", help = "sam file to convert")
parser.add_argument('-o', '--output_prefix', type = str,
                    nargs = "?", help = "a name to prefix the output files")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.sam_file is None or args.output_prefix is None:
    print("\n** required input missing\n"
          "** a sam file and output prefix are all required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

sam_fl = args.sam_file
output_prefix = args.output_prefix

# read file initially to grab paired / single names

print("Scanning sam file and building list of names...")

paired_set = set()
single_set = set()

with open(sam_fl) as fl:
    for line in fl:
        line = line.strip()
        if line.startswith("@"):
            continue
        name = line.split()[0]
        if name in single_set:
            paired_set.add(name)
            single_set.remove(name)
        else:
            single_set.add(name)

print("detected {} paired reads and {} single reads".format(len(paired_set), len(single_set)))

# only create files for required read types
# i.e. if there are no single end reads, no need to create empty file

if len(paired_set) > 0:
    forward_handle = open(args.output_prefix + ".R1.fastq", "w")
    reverse_handle = open(args.output_prefix + ".R2.fastq", "w")

if len(single_set) > 0:
    orphan_handle = open(args.output_prefix + ".single.fastq", "w")

print("writing output files")

with open(sam_fl) as fl:
    first_pair = True
    for line in fl:
        line = line.strip()
        if line.startswith("@"):
            continue
        cols = line.split("\t")
        name = cols[0]
        seq = cols[9]
        qual = cols[10]
        if name in paired_set:
            if first_pair == True:
                forward_handle.write("@{}/1\n{}\n+\n{}\n".format(name, seq, qual))
                first_pair = False
            else:
                reverse_handle.write("@{}/2\n{}\n+\n{}\n".format(name, seq, qual))
                first_pair = True
        else:
            orphan_handle.write("@{}\n{}\n+\n{}\n".format(name, seq, qual))















