#!/usr/bin/env python3

import argparse
from Bio import Phylo

parser = argparse.ArgumentParser(description="Given a tree(s) and a table counts the number of unique attributes and prints to the console.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", nargs='+', required=True,
						  help="A tree file(s).")

requiredArgs.add_argument("-a", "--attributes", dest="attributes", required=True,
						  help="A tab separated table containing two columns for the tip names and the given attribute (if more than 2 columns, the rest will be ignored).")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format. Accepted formats are: newick (default), nexus, nexml, phyloxml or cdao. If more than one tree, all trees should have the same format.")

args = parser.parse_args()

attributes = {}
count = {}
for line in open(args.attributes):
	tmp = line.strip().split()
	attributes[tmp[0]]=tmp[1]
	if tmp[1] not in count:
		count[tmp[1]] = 0

for filei in args.tree:
	# Read the phylogenetic tree
	T = Phylo.read(filei, args.formaTree)
	# Reset the count
	for k in count.keys():
		count[k] = 0
	# Count the number of attributes
	missed = 0
	for tip in T.get_terminals():
		if tip.name in attributes:
			tmp = attributes[tip.name]
			count[tmp] += 1
		else:
			missed += 0
	# print the number of attributes
	t = T.count_terminals()
	for a, c in count.items():
		p = round(c/t, 4)
		if len(args.tree) == 1:
			print(a, "\t", c, "\t", p, sep="")
		else:
			print(filei, "\t", a, "\t", c, "\t", p, sep="")
	# print if tips without attributes
	if missed != 0:
		if len(args.tree) == 1:
			print(missed, "tips had no given attribute")
		else:
			print(filei, "\t",missed, "tips had no given attribute")
