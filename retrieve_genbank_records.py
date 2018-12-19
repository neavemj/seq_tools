#!/usr/bin/env python

# reteive complete genbank records given a file of accession numbers
# Matthew J. Neave 19.12.2018

# requires biopython
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "matthewjneave1@gmail.com"

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("download complete genbank records given file containing NCBI accession numbers (one per line). "
                                 "Each record will be written into a separate file named by the accession number + .gb"

parser.add_argument('-l', '--acc_list', type = str,
        nargs = "?", help = "file containing list of NCBI accessions (one per line)")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.acc_list is None:
    print("\n** required input missing\n"
          "** a file containing wanted accession numbers is required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

acc_file = args.acc_list

# function to fetch and write out a genbank record

def retrieve_ncbi_record(ncbi_id):
    new_handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="genbank")
    seq_record = SeqIO.read(new_handle, "genbank")
    SeqIO.write(seq_record, ncbi_id + ".gb", "genbank")
