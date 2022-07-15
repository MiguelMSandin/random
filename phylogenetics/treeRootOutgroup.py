#!/usr/bin/env python3

import argparse
from ete3 import Tree

parser = argparse.ArgumentParser(description="Roots a newick tree by creating an outgropup on the last common ancestor of given tip names.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file in newick format.")

parser.add_argument("-o", "--output", dest="out", required=False, action="store",
					help="The output tree file. If not selected, will add '_rooted.tre' before the extension of the input tree file.")

parser.add_argument("-l", "--list", dest="list", required=False, default=None,
					help="A list of tip names where each tip is in one line.")

parser.add_argument("-n", "--names", dest="names", required=False, nargs='+', default=None,
					help="The name of the tips to be used for rooting, separated by a space. If selected, it will ignore the list given with the option '-l/--list'. Recommended for large trees.")

parser.add_argument("-p", "--pattern", dest="pattern", required=False, nargs='+', default=None,
					help="A text pattern(s) to be searched among the tip names. It will add all the tip names found to those given by either the 'names' or the 'list' options, if selected.")

parser.add_argument("-i", "--ignore", dest="ignore", required=False, action="store_true",
					help="If selected, will ignore the given tip names not present in the tree and only used those found.")

parser.add_argument("-s", "--sort", dest="sort", required=False, action="store_false",
					help="If selected, will not sort the tree branches according to the number of tips.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# setting variables --------------------------------------------------------------------------------
if args.out is None:
	import re
	out = re.sub("\\.[^\\.]+$", "_rooted.", args.tree) + re.sub(".*\\.", "", args.tree)
else:
	out = args.out

# Reading the files --------------------------------------------------------------------------------
if args.verbose:
	print("  Reading tree file:", args.tree)
T = Tree(args.tree, format=1)

# Reading the tips for rooting
tips = list()
if args.list is not None and args.names is None:
	tips = [line.strip() for line in open(args.list)]
if args.names is not None:
	if args.list is not None:
		print("  Option '-l/--list' is ignored")
	for tip in args.names:
		tips.append(tip)
if args.pattern is not None:
	for leaf in T.get_leaves():
		for pattern in args.pattern:
			if pattern in leaf.name:
				tips.append(leaf.name)

if args.names is None and args.list is None and args.pattern is None:
	import sys
	print("\nWarning! You have to select tips for rooting by either a file containing the tips with '-l/--list', by directly typing the tips with the option '-n/--names' or by providing a text pattern to be searched for\n")
	sys.exit(1)

# Ignoring tips if selected ------------------------------------------------------------------------
if args.ignore:
	count = 0
	for tip in tips:
		if tip not in T.iter_leaf_names():
			count += 1
			tips.remove(tip)
	if args.verbose:
		print("  In total", count, "tips were not found in the tree and therefore ignored")

# Rooting ------------------------------------------------------------------------------------------
if args.verbose:
	print("  Rooting tree", end="")
ancestor = T.get_common_ancestor(tips)
if args.verbose:
	print("\r  Rooting tree with a total of", len(ancestor), "tips")
T.set_outgroup(ancestor)

# Sorting ------------------------------------------------------------------------------------------
if args.sort:
	if args.verbose:
		print("  Sorting")
	T.ladderize()

# Writing file -------------------------------------------------------------------------------------
if args.verbose:
	print("  writing rooted tree file to:", out)
T.write(format=1, outfile=out)

if args.verbose:
	print("Done")
