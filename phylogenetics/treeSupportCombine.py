#!/usr/bin/env python3

import argparse
from Bio import Phylo

parser = argparse.ArgumentParser(description="Will combine the node support from different trees into the first input tree in the given order, and export an annotated tree.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--trees", dest="trees", required=True, nargs="+",
						  help="Tree files in newick format.")

# parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
# 					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")
# 
# parser.add_argument("-s", "--supportName", dest="supportName", required=False, default=None,
# 					help="The label given to the support value (e.g.: 'bootstrap', 'bs', 'posteriorp', ...)")

parser.add_argument("-d", "--delimiter", dest="delimiter", required=False, action="store", default="|",
					help="The delimiter to mark different confidence values. Default '|'.")

parser.add_argument("-a", "--absent", dest="absent", required=False, action="store", default="-",
					help="The character to represent non-existing clades. Default '-'.")

parser.add_argument("-o", "--output", dest="output", required=False, action="store", default=None,
					help="The output file name. By default will add '_condifenceCombined' to the first input tree with the given extension.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Setting output name and variables ----------------------------------------------------------------
if args.output is None:
	import re
	output = re.sub("\\.[^\\.]+$", "_condifenceCombined.", args.trees[0]) + re.sub(".*\\.", "", args.trees[0])
else:
	output = args.output

delim = str(args.delimiter)
absent = str(args.absent)

if len(args.trees) == 1:
	print(  "Warning, only one tree has been given.", flush=True)

# Reading files and checking internal nodes --------------------------------------------------------
i = 0
for tree in args.trees:
	i += 1
	if i == 1:
		if args.verbose:
			print("  Reading tree", tree, flush=True)
		# T = Phylo.read(tree, args.formaTree)
		T = Phylo.read(tree, "newick")
		nodesSupport = {}
		nodesTips = {}
		c = 0
		for clade in T.get_nonterminals():
			c += 1
			nodesSupport["node" + str(c)] = clade.confidence
			names = []
			for tip in clade.get_terminals():
				names.append(tip.name)
			nodesTips["node" + str(c)] = str(names)
	else:
		if args.verbose:
			print("\r  Reading tree ", i, "/", len(args.trees), " and referencing clades", sep="", end="", flush=True)
		# Ti = Phylo.read(tree, args.formaTree)
		Ti = Phylo.read(tree, "newick")
		nodesSupporti = {}
		nodesTipsi = {}
		c = 0
		for clade in Ti.get_nonterminals():
			c += 1
			nodesSupporti["node" + str(c)] = clade.confidence
			names = []
			for tip in clade.get_terminals():
				names.append(tip.name)
			nodesTipsi["node" + str(c)] = str(names)
		for node, tips in nodesTips.items():
			found = None
			for nodei, tipsi in nodesTipsi.items():
				if sorted(tips) == sorted(tipsi):
					found = nodei
			if found is not None:
				nodesSupport[node] = str(nodesSupport[node]) + str(delim) + str(nodesSupporti[found])
			else:
				nodesSupport[node] = str(nodesSupport[node]) + str(delim) + str(absent)

# Annotating output file ---------------------------------------------------------------------------
if args.verbose:
	print("\n  Annotating tree", flush=True)
c = 0
for clade in T.get_nonterminals():
	c += 1
	clade.comment = '&combined_confidence="' + nodesSupport["node" + str(c)] + '"'

# Writing file -------------------------------------------------------------------------------------
if args.verbose:
	print("  Exporting tree to", output, flush=True)
Phylo.write(T, output, "nexus")

if args.verbose:
	print("Done", flush=True)
