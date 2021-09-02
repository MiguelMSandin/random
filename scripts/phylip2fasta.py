#!/usr/bin/env python3

import argparse
import sys


parser = argparse.ArgumentParser(description="Converts a phylip file to fasta format, respecting the sequence name length.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input phylip file. Remember to replace spaces in sequence name beforehand, or the sequence name will be truncated at the spaces.")

requiredArgs.add_argument("-o", "--output", dest="file_out", required=True,
                    help="Output fasta file.")

args = parser.parse_args()

with open(args.file_out, "w") as outfile:
	j = 0
	for i in open(args.file_in):
		j += 1
		line = i.strip().split()
		name = line[0]
		sequence = line[1]
		if j == 1:
			alignLength = int(sequence)
			print("Alignment has", name, "sequences and ", sequence, "positions\nConverting")
		else:
			if len(sequence) != alignLength:
				print("\nERROR: Alignment contains sequences of different length")
				print("  Sequence", j, "has", len(sequence), "positions")
				print("Stop converting\n")
				sys.exit(1)
			else:
				print(">", name, "\n", sequence, sep="", file=outfile)

