#!/usr/bin/env python3

import argparse
from ete3 import Tree
import sys

parser = argparse.ArgumentParser(description="Computes the topological distance from a file containting 2 or more newick tree files and returns a pair-wise list of the two compared trees and the distance in a three columns tab separated table.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--trees", dest="trees", required=True,
						  help="A file containining 2 or more phylogenetic trees in newick format.")

parser.add_argument("-o", "--output", dest="out", required=False, action="store",
					help="The output file name. If not selected, will replace the extension of the input tree file by '_distances.list'.")

parser.add_argument("-n", "--names", dest="names", required=False, default=None,
					help="The list of tree names. Should be the same length as the given tree file.")

parser.add_argument("-a", "--all", dest="all", required=False, action="store_true",
					help="If selected, will compare all against all, without avoiding reciprocals (A-B and B-A) and self hits.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# setting variables --------------------------------------------------------------------------------
if args.out is None:
	import re
	out = re.sub("\\.[^\\.]+$", "_distances.list", args.trees)
else:
	out = args.out

# Reading the files --------------------------------------------------------------------------------
if args.verbose:
	print("  Reading tree file:", args.trees)
#T = Tree(args.trees, format=1)
tree_collection = []
for newick in open(args.trees):
	tree_collection.append(Tree(newick))

if args.names is not None:
	names = list()
	for line in open(args.names):
		names.append(line.split("\n")[0])
	if len(names) != len(tree_collection):
		print("  Warning! Number of given names (", len(names), ") do not match number of given trees (", len(tree_collection), ")", sep="")
		sys.exit(1)

# creating a list for the pair-wise comparison
pairs = set()
for i in range(0, len(tree_collection)):
	for j in range(0, len(tree_collection)):
		if args.all:
			if str(str(i) + "-" + str(j)) not in pairs:
				pairs.add(str(str(i) + "-" + str(j)))
		else:
			if i != j and str(str(i) + "-" + str(j)) not in pairs and str(str(j) + "-" + str(i)) not in pairs:
				pairs.add(str(str(i) + "-" + str(j)))

# Coparing all pair-wise pairs ---------------------------------------------------------------------
if args.verbose:
	i = 0
	print("  Comparing", str(len(tree_collection)), "trees in a total of", str(len(pairs)), "pair-wise comparisons")
with open(out, "w") as outFile:
	for pw in pairs:
		if args.verbose:
			i += 1
			print("\r    ", str(round(i / len(pairs) * 100)), "%", sep="", end="")
		t1 = int(pw.split('-')[0])
		t2 = int(pw.split('-')[1])
		rf=tree_collection[t1].robinson_foulds(tree_collection[t2], unrooted_trees=True)[0]
		if args.names is not None:
			print(str(names[t1]) + "\t" + str(names[t2]) + "\t" + str(rf), file=outFile)
		else:
			print(str(t1+1) + "\t" + str(t2+1) + "\t" + str(rf), file=outFile)

# Writing file -------------------------------------------------------------------------------------
if args.verbose:
	print("\n  List of pair-wise comparisons written to:", out)

if args.verbose:
	print("Done")
