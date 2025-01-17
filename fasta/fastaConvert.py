#!/usr/bin/env python3

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Changes the names of a fasta file given a tab delimited table with the old names in one column and the new names in the second column.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-i", "--input", dest="file_in", required=True,
                    help="Input fasta file.")

requiredArgs.add_argument("-o", "--output", dest="file_out", required=False,
                    help="Output fasta file.")

parser.add_argument("-f", "--from", dest="formatIn", required=False,
					default="fasta",
                    help="Input format. By default = %(default)s.")

parser.add_argument("-t", "--to", dest="formatOut", required=False,
					default="stockholm",
                    help="Output format. By default = %(default)s.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Reading table ------------------------------------------------------------------------------------
if args.verbose:
	print("  Converting '", args.file_in, "' from '", args.formatIn, "' to '", args.formatOut, "'", sep="")
SeqIO.convert(args.file_in, args.formatIn, args.file_out, args.formatOut)

if args.verbose:
	print("  Exported file:", args.file_out)
	print("Done")
