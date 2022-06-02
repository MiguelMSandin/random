#!/usr/bin/env python3

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Converts a fasta file where the sequences are in multiple lines to a fasta file where each sequence is in one line.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input fasta file.")

requiredArgs.add_argument("-o", "--output", dest="file_out", required=True,
                    help="Output fasta file.")

args = parser.parse_args()

with open(args.file_out, "w") as outfile:
    for line in SeqIO.parse(open(args.file_in), "fasta"):
        print(">" + line.id + "\n" + line.seq, file=outfile)

