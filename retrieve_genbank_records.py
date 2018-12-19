#!/usr/bin/env python

# reteive complete genbank records given a file of accession numbers
# Matthew J. Neave 19.12.2018

# requires biopython
from Bio import Entrez
from Bio import SeqIO
import argparse
import sys

Entrez.email = "matthewjneave1@gmail.com"

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("Given file containing NCBI accession numbers (one per line), download complete genbank records.\n"
                                 "Each record will be written into a separate file named by the accession number + .gb.\n"
                                 "To run directly without a file, do 'retrieve_ncbi_record.py -l <(echo MG924986)'\n")

parser.add_argument('-l', '--acc_file', type = str,
        nargs = "?", help = "file containing list of NCBI accessions (one per line)")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.acc_file is None:
    print("\n** required input missing\n"
          "** a file containing wanted accession numbers is required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

acc_handle = args.acc_file

# function to fetch and a genbank record

def retrieve_ncbi_record(ncbi_id):
    print("retrieving {} from NCBI".format(ncbi_id))
    new_handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="genbank")
    seq_record = SeqIO.read(new_handle, "genbank")
    return(seq_record)

# get list of accession numbers and Processing

acc_list = [line.strip() for line in open(acc_handle)]

for acc in acc_list:
    record = retrieve_ncbi_record(acc)
    SeqIO.write(record, acc + ".gb", "genbank")
