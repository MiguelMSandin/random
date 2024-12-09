#!/usr/bin/env python3

import argparse
import re

parser = argparse.ArgumentParser(description="Creates a network file from a matrix of similarities.O")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
					help="A tab separated matrix file, with headers and row names.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
					help="Output file. By default will add '_net.net' to the input file")

parser.add_argument("-t", "--threshold", dest="threshold", required=False, default=1, type=float,
					help="Considers only hits equal or above the given threshold. By default=1")

parser.add_argument("-a", "--addHeaders", dest="addHeaders", required=False, action="store_true",
					help="If selected, will add the following headers to the net: 'source target id'.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

if args.file_out is None:
	out = re.sub("\\.[^\\.]+$", "_net.net", args.file_in)
else:
	out = args.file_out

# Reading file
if args.verbose:
	print("  Reading file")
with open(out, "w") as outfile:
	if args.addHeaders:
		print("source\ttarget\tid", file=outfile)
	names = list()
	i = 0
	for line in open(args.file_in):
		i += 1
		linei = line.strip().split('\t')
		if i == 1:
			j = 0
			for tmp in linei:
				j += 1
				if j != 1:
					names.append(tmp)
		else:
			hit = linei[0]
			hits = linei[1:]
			hitc = 0
			for h in hits:
				hitc += 1
				if float(h) >= args.threshold:
					tmp = names[hitc-1]
					print(hit, "\t", tmp, "\t", h, sep="", file=outfile)

if args.verbose:
	print("Done")
