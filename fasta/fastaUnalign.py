#!/usr/bin/env python3

import argparse
import re

parser = argparse.ArgumentParser(description="Removes gaps from sequences in a fasta file.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="fileIn", nargs='+', required=True,
						  help="Input fasta file(s).")

parser.add_argument("-o", "--output", dest="fileOut", nargs='+', required=False, 
					default=None,
					help="Output fasta file(s). By default will append '_unaligned' to the end of the file name (respecting the extension).")

parser.add_argument("-g", "--gaps", dest="gap", required=False, 
					default="-",
					help="The character representing the gaps. By default='-'.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

if args.fileOut is not None:
	i = 0
	if len(args.fileIn) != len(args.fileOut):
		import sys
		print("Error: Number of input files do not match number of output files.")
		sys.exit(1)

for filei in args.fileIn:
	if args.verbose:
		print("  Unaligning", filei)
	if args.fileOut is None:
		outFile = re.sub("\\.[^\\.]+$", "_unaligned.", filei) + re.sub(".*\\.", "", filei)
	else:
		outFile = args.fileOut[i]
		i += 1
	with open(outFile, "w") as outfile:
		for line in open(filei):
			if line.startswith(">"):
				print(line, end="", file=outfile)
			else:
				lineout = re.sub(args.gap, "", line)
				print(lineout, end="", file=outfile)

if args.verbose:
	print("Done")
