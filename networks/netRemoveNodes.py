#!/usr/bin/env python3

import argparse
import re

parser = argparse.ArgumentParser(description="Removes nodes from a network either with a given pattern or selected in a list.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True, nargs="+",
                    help="Input file. A tab separated file containing at least 2 columns (an edge list), the rest of the columns are irrelevant.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None, nargs="+",
                    help="Output file. By default will add '_clean' to the input file name (respecting the extension).")

parser.add_argument("-l", "--list", dest="listn", required=False, default=None,
                    help="A list of nodes, with each node in one line.")

parser.add_argument("-p", "--pattern", dest="pattern", required=False, default=None, nargs="+",
                    help="A text string contained in the nodes to look for.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

if (args.listn is None) & (args.pattern is None):
	print("  Error: Please select the nodes to be remove with with a list and/or a pettern.\nExiting")
	import sys
	sys.exit(1)

if args. file_out is not None:
	if len(args.file_in) != len(args.file_out):
		print("  Error: The number of input files (", len(args.file_in), ") is different than the number of output files (", len(args.file_out), ") provided\nExiting", sep="")
		import sys
		sys.exit(1)

toRemove = set()
if args.listn is not None:
	for line in open(args.listn):
		line = re.sub("\n$", "", line)
		toRemove.add(line)
	if args.verbose:
		print("  In total,", len(toRemove), "nodes will be removed")

i = 0
for filei in args.file_in:
	# Output file
	if args.file_out is None:
		outFile = re.sub("\\.[^\\.]+$", "_clean", args.file_in[i]) + re.sub(".*\\.", ".", args.file_in[i])
	else:
		outFile = args.file_out[i]
	if args.verbose:
		print("  Reading network:", filei)
	edgesin = 0
	edgesout = 0
	with open(outFile, "w") as outfile:
		for line in open(filei):
			export = True
			edgesin += 1
			linei = line.strip().split()
			node1 = linei[0]
			node2 = linei[1]
			if args.listn is not None:
				if (node1 in toRemove) or (node2 in toRemove):
					export = False
			if args.pattern is not None:
				for pattern in args.pattern:
					if (pattern in node1) or (pattern in node2):
						export = False
			if export:
				edgesout += 1
				lineout = re.sub("\n$", "", line)
				print(lineout, file=outfile)
		if args.verbose:
			print("    Edges in:  ", edgesin)
			print("    Edges out: ", edgesout)
			print("    Clean network exported to:", outFile)
	i += 1
