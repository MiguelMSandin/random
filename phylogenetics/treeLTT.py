#!/usr/bin/env python3

import argparse
#from Bio import Phylo
import treeswift

parser = argparse.ArgumentParser(description="From a time calibrated tree, it will export a tab delimited table with the number of Lineages Through Time (LTT).")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file in newick format.")

parser.add_argument("-o", "--out", dest="out", required=False, default=None,
						  help="The output file name. By default, will remove the extension of the input tree file and add '_LTT.tsv'.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default) and nexus.")

parser.add_argument("-s", "--subtrees", dest="subtrees", required=False, action="store_true",
					help="If selected, will also append to the table the LTT data for every annotated clade. NOT IMPLEMENTED YET!")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Setting output name tree -------------------------------------------------------------------------
if args.out is None:
	import re
	out = re.sub("\\.[^\\.]+$", "_LTT.tsv", args.tree)
else:
	out = args.out

# Reading file -------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading tree file:", args.tree)
if args.formaTree == "newick":
	tree = treeswift.read_tree_newick(args.tree)
elif args.formaTree == "nexus":
	tree = treeswift.read_tree_nexus(args.tree)
	tree = tree['tree_1']
else:
	import sys
	print("Error: Accepted formats are 'newick' or 'nexis'.")
	sys.exit(1)

# Extracting LTT -----------------------------------------------------------------------------------
if args.verbose:
	print("  Calculating Lineages Through Time")
table = tree.lineages_through_time(show_plot=False)

# Exporting ----------------------------------------------------------------------------------------
if args.verbose:
	print("  Exporting table to:", out)
with open(out, "w") as outfile:
	print("time\tlineages", file=outfile)
	for time, lineages in table.items():
		print(str(time) + "\t" + str(lineages), file=outfile)

if args.verbose:
	print("Done")
