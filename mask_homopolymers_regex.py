#!/usr/bin/env python
# python 3

# script to mask homopolymers (remove) in a multiple sequence alignment
# Matthew J. Neave 29.09.2020 <matthewjneave1@gmail.com>

# library imports

import sys
import argparse
import re
from Bio import AlignIO # requires biopython

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("Mask homopolymer regions (delete) in a multiple sequence alignment")

parser.add_argument('-a', '--alignment', type = str,
        nargs = "?", help = "alignment in fasta format")
parser.add_argument('-l', '--homopolymer_length', type = str,
        nargs = "?", help = "length of homopolymers to be removed")
parser.add_argument('-f', '--flanking_length', type = int,
        nargs = "?", help = "length of flanking sequence to be removed")
parser.add_argument('-o', '--output', type = str,
        nargs = "?", help = "output file name")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.alignment is None or \
    args.homopolymer_length is None or \
    args.flanking_length is None or \
    args.output is None:
    print("\n** required input missing\n"
          "** an alignment file, minimum homopolymer length and output name are all required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)


alignment_file = args.alignment
output_file = args.output
homopolymer_length = args.homopolymer_length
flanking_length = args.flanking_length

# iterate through alignment positions

align = AlignIO.read(alignment_file, "fasta")
cuts = []

for nucl in ["A", "C", "T", "G"]:
    # regex to match nucleotide of homopolymer_length or greater
    pat = re.compile(nucl + "{" + homopolymer_length + ",}")
    # just search for homopolymers in the first sequence
    # add flanking nucs and trim
    # could cause something weird to happen but this is easiest
    # could consider adding option to indicate ref seq and use this for homopolymer identification
    ref_seq = align[0]
    for match in pat.finditer(str(ref_seq.seq).upper()):
        # add flanking length to cut sites
        start_cut = match.start() - flanking_length
        end_cut = match.end() + flanking_length
        cuts.append([start_cut, end_cut]) 


# need to sort cuts so that I can build up the new alignment bit by bit
cuts_sorted = sorted(cuts)
align_cuts = []
count = 0

# need enumerate because I go back to previous list item for lower cut site
for index, cut in enumerate(cuts_sorted):
    count += 1
    print("** Masking homopolymer region: {} - {} **\n".format(cut[0], cut[1]), align[:, cut[0]:cut[1]])
    if count == 1: 
        # first alignment bit is from 0 to the first start cut site
        #print("cutting: ", 0, cut[0])
        masked_align = align[:, 0:cut[0]]
    elif count == len(cuts_sorted):
        # final alignment bit is from previous end site to this cut site, PLUS
        #print("cutting: ", cuts_sorted[index - 1][1],cut[0])
        masked_align = masked_align + align[:, cuts_sorted[index - 1][1]:cut[0]]
        # final bit is last cut site to end of alignment
        #print("cutting: ", cut[1], "end")
        masked_align = masked_align + align[:, cut[1]:]
    else:
        # middle alignment bits are from the previous end site to the current start site
        #print("cutting: ", cuts_sorted[index - 1][1],cut[0])
        masked_align = masked_align + align[:, cuts_sorted[index - 1][1]:cut[0]]
    #print(masked_align.get_alignment_length())

# now write newly built alignment
AlignIO.write(masked_align, output_file, "fasta")

# print some output

print("\n** Removed {} homopolymer regions from {} and saved to {}".format(count, alignment_file, output_file))


