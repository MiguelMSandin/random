#!/usr/bin/env python3

import argparse
import re

parser = argparse.ArgumentParser(description="Removes reciprocal hits (A-B = B-A) and self (A=A) hits from a tsv table where the first two columns are the identifiers.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input file. This assumes a file with the following columns: 'qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen'.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
                    help="Output file. Returns the filtered file to the specified location. By default will add '_clean.net' to the input file, after removing the extension.")

parser.add_argument("-r", "--reciprocal", dest="reciprocal", required=False, action="store_true",
					help="If selected, will not remove reciprocal hits.")

parser.add_argument("-s", "--self", dest="selfHits", required=False, action="store_true",
					help="If selected, will not remove self hits.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

if args.file_out is None:
	out = re.sub("\\.[^\\.]+$", "_clean.net", args.file_in)
else:
	out = args.out

if args.reciprocal:
	endr = "(not removed)"
	if args.verbose:
		print("Ignoring reciprocal hits")
else:
	endr = ""

if args.selfHits:
	ends = "(not removed)"
	if args.verbose:
		print("Ignoring self hits")
else:
	ends = ""

if args.verbose:
	print("Cleaning")

clean = set()
lines = 0
reciprocal = 0
selfHit = 0
accepted = 0
with open(out, "w") as outfile:
	for line in open(args.file_in):
		toPrint = True
		lines += 1
		linei = line.strip().split()
		seq1 = linei[0]
		seq2 = linei[1]
		hit = str(seq1) + "-" + str(seq2)
		hitr =str(seq2) + "-" + str(seq1)
		if seq1 == seq2:
			selfHit += 1
			toPrint = False
			if args.selfHits:
				toPrint = True
				accepted += 1
		elif hit in clean or hitr in clean:
			reciprocal += 1
			toPrint = False
			if args.reciprocal:
				toPrint = True
				accepted += 1
		else:
			accepted += 1
			clean.add(str(hit))
		if toPrint:
			print(line, file=outfile, end="")

if args.verbose:
	print("  Hits in:         ", lines)
	print("  Reciprocal hits: ", reciprocal, endr)
	print("  Self hits:       ", selfHit, ends)
	print("  Hits out:        ", accepted)
	print("Done")
