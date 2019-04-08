#!/usr/bin/env python

### extract gene name and protein sequence from genbank files ##
# Matthew J Neave. 26.9.14 #

import sys
from Bio import SeqIO


gbk_filename = sys.argv[1]
faa_filename = sys.argv[2]
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print "Analysing GenBank record %s" % seq_record.id
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            output_handle.write(">%s %s\n%s\n" % (
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_feature.qualifiers['product'][0],
                   seq_feature.qualifiers['translation'][0]))

output_handle.close()
input_handle.close()
