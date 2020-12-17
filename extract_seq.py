#!/usr/bin/python3

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Extracting sequences from a fasta file (-f) that match a pattern (-p) and exporting to output file (-o)")

# Add the arguments to the parser
parser.add_argument("-p", "--patern", dest="patern", required=True,
                    help="pattern to be matched (regex)")
parser.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input file")
parser.add_argument("-o", "--output", dest="file_out", required=True,
                    help="Output file")
args = parser.parse_args()

with open(args.file_out, "w") as outfile:
    for i in SeqIO.parse(open(args.file_in), "fasta"):
        if args.patern in i.id:
            SeqIO.write([i], outfile, "fasta")
