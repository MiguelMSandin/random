#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import re

parser = argparse.ArgumentParser(description="Remove sequences from a fasta file (-f) that listed in a list (-l) and exporting to output file (-o)")

# Add the arguments to the parser
parser.add_argument("-l", "--list", dest="listSeq", required=True,
                    help="pattern to be matched (regex)")
parser.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input file")
parser.add_argument("-o", "--output", dest="file_out", required=True,
                    help="Output file")
args = parser.parse_args()

toRemove = [line.strip() for line in open(args.listSeq)]

file_outr = re.sub("\..*$", "_removed.fasta", args.file_out)

for i in SeqIO.parse(open(args.file_in), "fasta"):
    if i.id not in toRemove:
        with open(args.file_out, "a") as outfile:
            SeqIO.write([i], outfile, "fasta")
    elif i.id in toRemove:
        with open(file_outr, "a") as outfiler:
            SeqIO.write([i], outfiler, "fasta")
    else:
        print("Sequence '", i , "' not found in the fasta file", sep="")


#with open(args.file_out, "w") as outfile:
#    for i in SeqIO.parse(open(args.file_in), "fasta"):
#        if i.id not in toRemove:
#            SeqIO.write([i], outfile, "fasta")  
#
#with open(file_outr, "w") as outfiler:
#    for i in SeqIO.parse(open(args.file_in), "fasta"):
#        if i.id in toRemove:
#            SeqIO.write([i], outfiler, "fasta")
