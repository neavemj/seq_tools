#!/usr/bin/env python
# python 3

# script to 're-pair' mates after trimming
# Matthew J. Neave 5.5.2015 <matthewjneave1@gmail.com>
# updated to py3 26.06.2018

# library imports

import sys
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator # requires Biopython

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("re-pairing mate pair sequences")

parser.add_argument('forward_reads', type = argparse.FileType("r"),
        nargs = "?", help = "fastq file containing forward R1 reads")
parser.add_argument('reverse_reads', type = str,
        nargs = "?", help = "fastq file containing reverse R2 reads")
parser.add_argument('output_prefix', type = str,
        nargs = "?", help = "a name to prefix the output files")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

forward_handle = open(args.output_prefix + ".R1.fq", "w")
reverse_handle = open(args.output_prefix + ".R2.fq", "w")
orphan_handle = open(args.output_prefix + ".orphan.fq", "w")

print("Scanning reverse file and building list of names...")
reverse_ids = set()
paired_ids = set()
for title, seq, qual in FastqGeneralIterator(open(args.reverse_reads)):
    reverse_ids.add(title.split()[0])

print("Processing forward file")

for title, seq, qual in FastqGeneralIterator(args.forward_reads):
    name = title.split()[0]
    if name in reverse_ids:
        # paired reads
        paired_ids.add(name)
        reverse_ids.remove(name) # saves a little memory
        forward_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    else:
        # orphan reads
        orphan_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
forward_handle.close()

print("Processing reverse file")

for title, seq, qual in FastqGeneralIterator(open(args.reverse_reads)):
    name = title.split()[0]
    if name in paired_ids:
        # paired reads
        reverse_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    else:
        # orphaned reads
        orphan_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
orphan_handle.close()
reverse_handle.close()
print("done")