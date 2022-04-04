#!/usr/bin/env python3

import argparse
from Bio import Phylo

parser = argparse.ArgumentParser(description="Prunes a tree from tips present in a list.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

parser.add_argument("-o", "--output", dest="out", required=False, action="store",
					help="The output tree file. If not selected, will add '_rooted.tre' before the extension of the input tree file.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-F", "--formatOut", dest="formaTreeOut", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-l", "--list", dest="list", required=False, default=None,
					help="A list of tip names where each tip is in one line.")

parser.add_argument("-n", "--names", dest="names", required=False, nargs='+', default=None,
					help="The name of the tips to be used for rooting, separated by a space. If selected, it will ignore the list given with the option '-l/--list'. Recommended for large trees.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# setting variables --------------------------------------------------------------------------------
if args.out is None:
	import re
	out = re.sub("\\.[^\\.]+$", "_rooted.", args.tree) + re.sub(".*\\.", "", args.tree)
else:
	out = args.out

# Reading the tips for rooting
if args.list is not None and args.names is None:
	tips = list()
	tips = [line.strip() for line in open(args.list)]
if args.names is not None:
	if args.list is not None:
		print("  Option '-l/--list' is ignored")
	tips = list()
	for tip in args.names:
		tips.append(tip)
if args.names is None and args.list is None:
	import sys
	print("\nWarning! You have to select tips for rooting by either a file containing the tips with '-l/--list' or directly typing the tips with the option '-n/--names'.\n")
	sys.exit(1)

# Reading the file ---------------------------------------------------------------------------------
if args.verbose:
	print("  Reading tree file:", args.tree)
T = Phylo.read(args.tree, args.formaTree)

# Rooting ------------------------------------------------------------------------------------------
if args.verbose:
	if args.list is not None and args.names is None:
		print("  Rooting tree with tips in list:", args.list)
	if args.names is not None:
		print("  Rooting tree with the follwing tips:")
		for tip in tips:
				print("  -", tip)
T.root_with_outgroup(tips)

# Writing file -------------------------------------------------------------------------------------
if args.verbose:
	print("  writing rooted tree file to:", out)
Phylo.write(T, out, args.formaTreeOut)

if args.verbose:
	print("Done")
