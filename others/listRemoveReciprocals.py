#!/usr/bin/env python3

import argparse
import re

parser = argparse.ArgumentParser(description="Removes reciprocal (A-B = B-A) and identical (A-A) hits from a tab separated list. The first and second columns are taken as identifiers. The rest of the columns are irrelevant.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input file.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
                    help="Output file. By default will add '_clean' to the input file name (respecting the extension).")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Output file
if args.file_out is None:
	outFile = re.sub("\\.[^\\.]+$", "_clean", args.file_in) + re.sub(".*\\.", ".", args.file_in)
else:
	outFile = args.file_out

if args.verbose:
	print("  Setting hit names")
clean = {}
entries = 0
reciprocal = 0
selfHit = 0
accepted = 0
for line in open(args.file_in):
	entries += 1
	seq = line.strip().split()
	seq1 = seq[0]
	seq2 = seq[1]	
	if seq1 == seq2:
		selfHit = selfHit + 1	
	hit = str(seq1) + "\t" + str(seq2)
	hitr =str(seq2) + "\t" + str(seq1)	
	if hit in clean:
		reciprocal += 1
	elif hitr in clean:
		reciprocal += 1
	else:
		clean[hit] = line
		accepted += 1

if args.verbose:
	print("  Writing cleaned file to: '", outFile, "'", sep="")
with open(outFile, "w") as outfile:
    for hit in list(clean.keys()):
        print(clean[hit], end="", file=outfile)

if args.verbose:
	print("    Entries:         ", entries, sep="")
	print("    Removed:         ", reciprocal+selfHit, " (", round((reciprocal+selfHit)/entries*100,2), "%)", sep="")
	print("      Reciprocal:    ", reciprocal, " (", round((reciprocal)/entries*100,2), "%)", sep="")
	print("      Self matching: ", selfHit, " (", round((selfHit)/entries*100,2), "%)", sep="")
	print("    Out:             ", accepted, " (", round(accepted/entries*100,2), "%)", sep="")
	print("Done")
