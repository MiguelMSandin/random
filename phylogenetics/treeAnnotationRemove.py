#!/usr/bin/env python3

import argparse
from Bio import Phylo
import re

parser = argparse.ArgumentParser(description="Given an annotated tree, removes all kind of annotations at nodes.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The annotated tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-o", "--output", dest="output", required=False, default=None,
					help="The output name of the tree file. By default will add '_clean' to the input file name and the original extension.")

parser.add_argument("-F", "--formatOut", dest="formaTreeOut", required=False, default='newick',
					help="The output tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

methodsHelp = parser.add_argument_group('Methods')

args = parser.parse_args()

# Reading files ------------------------------------------------------------------------------------
if args.verbose:
	print("Reading annotated tree:  ", args.tree)
T = Phylo.read(args.tree, args.formaTree)

# Setting variables --------------------------------------------------------------------------------

# Output file
if args.output is None:
	outFile = re.sub("\\.[^\\.]+$", "_clean.", args.tree) + re.sub(".*\\.", "", args.tree)
else:
	outFile = args.output

# Loop through internal nodes and look for childs -------------------------------------------------
if args.verbose:
	print("  Removing annotations")

for line in T.get_nonterminals():
	if line.name is not None:
		line.name = None
	elif line.comment is not None:
		line.comment = None


# Writing files ------------------------------------------------------------------------------------
if args.verbose:
	print("  Exporting cleaned tree to:", outFile)
Phylo.write(T, outFile, args.formaTreeOut)

if args.verbose:
	print("Done")
