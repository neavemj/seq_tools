#!/usr/bin/env python

### extract nucleotide or protein sequence from genbank files ##
# Matthew J Neave. 08.04.19 #

import sys
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser("Given NCBI genbank record return fasta files of nucl and prot seqs.\n"
                                 "Each sequence will be written into a file named by the accession number + .faa or .fna\n"
                                 "Only provide an output filename for the type you want\n"
                                 "For example, if you only provide a protein output, the nucleotides won't be written")

parser.add_argument('-g', '--genbank_file', type = str,
                    nargs = "?", help = "genbankfile")
parser.add_argument('-n', '--nucl_output', type = str,
                    nargs = "?", help = "fastq file containing reverse R2 reads")
parser.add_argument('-p', '--prot_output', type = str,
                    nargs = "?", help = "fastq file containing reverse R2 reads")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.genbank_file is None:
    print("\n** required input missing\n"
          "** a genbank file is required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)


for seq_record in SeqIO.parse(args.genbank_file, "genbank"):
    print("Analysing GenBank record {}".format(seq_record.id))
    if args.nucl_output:
        nucl_handle = open(args.nucl_output, "w")
        nucl_handle.write(">{} {}\n{}\n".format(
                seq_record.id,
                seq_record.description,
                seq_record.seq))

    if args.prot_output:
        prot_handle = open(args.prot_output, "w")
        for seq_feature in seq_record.features:
            if seq_feature.type=="CDS" :
                assert len(seq_feature.qualifiers['translation'])==1
                prot_handle.write(">{} {}\n{}\n".format(
                       seq_feature.qualifiers['protein_id'][0],
                       seq_feature.qualifiers['product'][0],
                       seq_feature.qualifiers['translation'][0]))

nucl_handle.close()
prot_handle.close()
