#!/usr/bin/env python3

import argparse
from Bio import Phylo

parser = argparse.ArgumentParser(description="Given a list of trees and a list of tip names, counts the number of tips in the Last Common Ancestor clade")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--trees", dest="trees", nargs='+', required=True,
						  help="A list of trees.")

requiredArgs.add_argument("-n", "--names", dest="names", nargs='+', required=True,
						  help="A character string(s) to be searched in the tree tip names separated with a space.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-s", "--summary", dest="summary", required=False, action="store_false",
					help="If selected, will not print a summary of the tip count at the end.")

parser.add_argument("-l", "--list", dest="list", required=False, default=None,
					help="If selected, will print to the given file all unique tip names found in the last common ancestor tip names.")

args = parser.parse_args()

n = list()
if args.list is not None:
	tips = {}
for tree in args.trees:
	T = Phylo.read(tree, args.formaTree)
	clade = T.common_ancestor(args.names)
	count = clade.count_terminals()
	print(str(count), "\t", tree)
	n.append(count)
	if args.list is not None:
		for tip in clade.get_terminals():
			if tip.name in tips:
				tips[tip.name] += 1
			else:
				tips[tip.name] = 1

if args.summary:
	import statistics as st
	import numpy as np
	print("")
	print("min\tmedian\tmax\t\tmean\tsd")
	mi = str(min(n))
	med = str(int(np.percentile(n, 50)))
	ma = str(max(n))
	mea = str(round(st.mean(n), 2))
	sd = str(round(st.stdev(n), 2))
	print(mi + "\t" + med + "\t" + ma + "\t\t" + mea + "\t" + sd)

if args.list is not None:
	print("\nWriting list of unique tips to:", args.list)
	with open(args.list, "w") as outlist:
		print("count\tselected\ttip", file=outlist)
		for tip, c in tips.items():
			if tip in args.names:
				s = "yes"
			else:
				s = "no"
			print(str(c) + "\t" + s + "\t" + str(tip), file=outlist)
	print("Done")
