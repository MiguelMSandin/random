#!/usr/bin/env python3

import argparse
from Bio import Phylo
import re

parser = argparse.ArgumentParser(description="From an ultrametric tree, extracts all node lengths (from node to the first tip of the given node) to a table.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file containing one or several trees. If several, all trees should have the same bifurcation patterns.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The annotated tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-o", "--output", dest="output", required=False, default=None, nargs="+",
					help="The output table name, with columns for the nodes and rows for the different heights(if more than one tree). By default will take the name of the input tree file name followed by '_nodeHeights.tsv'.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

methodsHelp = parser.add_argument_group('Methods')

args = parser.parse_args()

# Setting options ----------------------------------------------------------------------------------
if args.output is None:
	outFile = re.sub("\\.[^\\.]+$", "_nodeHeights.tsv", args.tree)
else:
	outFile = args.output

# Reading files ------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading input tree: ", args.tree)
trees = Phylo.parse(args.tree, args.formaTree)

with open(args.tree, 'r') as tc:
    for count, line in enumerate(tc):
        pass
treeCount = count + 1

# Exporting heights to table -----------------------------------------------------------------------
if args.verbose:
	print("  Exporting table to: ", outFile, end="")
i = 0
with open(outFile, 'w') as outfile:
	for tree in trees:
		heights = []
		i += 1
		print("\r  Exporting table to:  ", outFile, "\t", str(round(i/treeCount*100)), "%", sep="", end="")
		if i == 1:
			nnodes = len(tree.get_nonterminals())
			cols = []
			for j in range(0,nnodes):
				cols.append("n"+str(j+1))
			row = '\t'.join(cols)
			print(str(row), file=outfile)
		for node in tree.get_nonterminals():
			height = None
			for tip in node.get_terminals():
				if height is None:
					height = node.distance(node, tip)
					heights.append(str(height))
		row = '\t'.join(heights)
		print(str(row), file=outfile)

if args.verbose:
	print("\nDone")
