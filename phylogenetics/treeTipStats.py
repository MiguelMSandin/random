#!/usr/bin/env python3

import argparse
from Bio import Phylo

parser = argparse.ArgumentParser(description="For every tip of a phylogenetic tree, exports a tab delimited table with the tip name, the branch length, the distance to the root and the number of nodes to the root.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, 
					default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-o", "--output", dest="output", required=False,
					default=None,
					help="The name of the output file. By default, will remove the extension and add '_tipStats.tsv'.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Setting output name ------------------------------------------------------------------------------
if args.output is None:
	import re
	output = re.sub("\\.[^\\.]+$", "_tipStats.tsv", args.tree)
else:
	output = args.output

# Reading file -------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading file:", args.tree)
T = Phylo.read(args.tree, args.formaTree)

# Number of tips
if args.verbose:
	ntips = T.count_terminals()
	i = 0
	print("  Getting stats", end="")
with open(output, "w", buffering=1) as outfile:
	print("tipName\tbranchLength\tdistRoot\tnodesRoot", file=outfile)
	for tip in T.get_terminals():
		if args.verbose:
			i += 1
			print("\r  Getting stats\t", round((i/ntips)*100), "%", sep="", end="")
		height = round(T.distance(tip), 8)
		nodes = len(T.get_path(tip))
		print(str(tip.name) + "\t" + str(tip.branch_length) + "\t" + str(height) + "\t" + str(nodes), file=outfile)

if args.verbose:
	print("\n  Table exported to:", output)
	print("Done")
