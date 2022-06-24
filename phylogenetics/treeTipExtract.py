#!/usr/bin/env python3

import argparse
from Bio import Phylo
import re

parser = argparse.ArgumentParser(description="")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-o", "--output", dest="output", required=False, default=None,
					help="The output file name. By default will add to the file name '_tipNames.tsv'")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Output file
if args.output is None:
	outFile = re.sub("\\.[^\\.]+$", "_tipNames.tsv", args.tree)
else:
	outFile = args.output

# Reading file -------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading tree file:", args.tree)
T = Phylo.read(args.tree, args.formaTree)

# Rename tree --------------------------------------------------------------------------------------
if args.verbose:
	print("  Extracting")
with open(outFile, 'w') as outfile:
	for line in T.get_terminals():
		print(str(line.name), file=outfile)

if args.verbose:
	print("Done")
