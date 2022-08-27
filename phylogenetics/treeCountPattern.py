#!/usr/bin/env python3

import argparse
from Bio import Phylo

parser = argparse.ArgumentParser(description="Given a tree, counts the number of tips and prints to the console.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

requiredArgs.add_argument("-p", "--pattern", dest="pattern", nargs='+', required=True,
						  help="A character string(s) to be searched in the tree tip names, if more than one separated with a space. You can link several patterns as one with the '+' character (e.g.; A B C+D E: will search for A, then B, then C and D and lastly E).")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-s", "--search", dest="search", required=False, default='tips',
					choices=['tips', 'clades'],
					help="Whether to search among the tree tips ('tips', default) or the number of tips in the last common ancestor clade gathering all tips with the given pattern ('clades').")

args = parser.parse_args()

T = Phylo.read(args.tree, args.formaTree)

for pattern in args.pattern:
	if args.search == "tips":
		count = 0
		if "+" in pattern:
			for p in pattern.split("+"):
				for tip in T.get_terminals():
					if p in tip.name:
						count += 1
		else:
			for tip in T.get_terminals():
				if pattern in tip.name:
					count += 1
	if args.search == "clades":
		tips = list()
		if "+" in pattern:
			for p in pattern.split("+"):
				for tip in T.get_terminals():
					if p in tip.name:
						tips.append(tip.name)
		else:
			for tip in T.get_terminals():
				if pattern in tip.name:
					tips.append(tip.name)
		if len(tips) > 0:
			clade = T.common_ancestor(tips[0], tips[len(tips)-1])
			count = clade.count_terminals()
		else:
			count = 0
	print(pattern, "\t", str(count))
