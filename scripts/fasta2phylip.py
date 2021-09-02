#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import re
import sys


parser = argparse.ArgumentParser(description="Converts a fasta file to phylip format, respecting the sequence name length")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input fasta file. Remember to replace spaces in sequence name beforehand, or the sequence name will be truncated at the spaces.")

requiredArgs.add_argument("-o", "--output", dest="file_out", required=True,
                    help="Output fasta file.")

args = parser.parse_args()

seqs = 0
for line in open(args.file_in):
	if ">" in line:
		seqs += 1

with open(args.file_out, "w") as outfile:
	j = 0
	for i in SeqIO.parse(open(args.file_in), "fasta"):
		j += 1
		sequence = i.seq
		name = i.id
		if j == 1:
			alignLength = len(sequence)
		if j == 1:
			print(seqs, alignLength, file=outfile)
			print(name, "  ", sequence, file=outfile)
		else:
			if len(sequence) != alignLength:
				print("\nERROR: Alignment contains sequences of different length")
				print("  Sequence", j, "has", len(sequence), "positions")
				print("Stop converting\n")
				sys.exit(1)
			else:
				print(name, "  ", sequence, file=outfile)
