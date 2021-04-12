#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Converts the cases of a fasta file from lower to upper or the other way around.")

parser.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input fasta file.")

parser.add_argument("-o", "--output", dest="file_out", required=True,
                    help="Output fasta file.")

parser.add_argument("-c", "--case", dest="case", required=True, choices=['u', 'upper', 'l', 'lower'],
                    help="Changes from lower to upper cases ('u' or 'upper') or the other way around ('l' or 'lower').")

args = parser.parse_args()

with open(args.file_out, "w") as outfile:
	for line in open(args.file_in):
		if ">" in line:
			print(line, end="", file=outfile)
		else:
			if args.case == "u" or args.case == "upper":
				print(line.upper(), end="", file=outfile)
			if args.case == "l" or args.case == "lower":
				print(line.lower(), end="", file=outfile)

