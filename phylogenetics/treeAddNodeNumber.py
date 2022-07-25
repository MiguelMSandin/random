#!/usr/bin/env python3

import argparse
from Bio import Phylo
import subprocess
import re

parser = argparse.ArgumentParser(description="Takes a tree and exports the same tree with the node numbers as comments.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
                    help="A tree file.")

parser.add_argument("-o", "--output", dest="out", required=False, action="store",
					help="The output tree file. If not selected, will add '_nodeNumbers.tre' before the extension of the input tree file.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-n", "--numbers", dest="numbers", required=False,
					choices=['order', 'tips', 'invert'], default='order',
					help="How to number the nodes: 'order' (default): 1 to N; 'tips': T tp T+N; 'invert': N to 1. Being 'N' the number of internal nodes and 'T' the number of tips")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# setting variables --------------------------------------------------------------------------------
if args.out is None:
	out = re.sub("\\.[^\\.]+$", "_nodeNumbers.tre", args.tree)
else:
	out = args.out

# Reading files ------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading tree")
T = Phylo.read(args.tree, args.formaTree)

numbers = list()
if args.numbers == "order":
	i = 0
	for t in T.get_nonterminals():
		i += 1
		numbers.append(i)
if args.numbers == "tips":
	i = T.count_terminals()
	for t in T.get_nonterminals():
		i += 1
		numbers.append(i)
if args.numbers == "invert":
	i = 0
	for t in T.get_nonterminals():
		i += 1
	numbers = list(reversed(range(i+1)))

# Adding node heights to node comments -------------------------------------------------------------
if args.verbose:
	print("  Annotating node numbers")
i = 0
for clade in T.get_nonterminals():
	if clade.comment is None:
		clade.comment = str("[&nodeNumber=" + str(numbers[i]) + "]")
	else:
		tmp = re.sub("\\]", "", clade.comment)
		clade.comment = str(tmp + ",nodeNumber=" + str(numbers[i]) + "]")
	i += 1

# Writing file -------------------------------------------------------------------------------------
if args.verbose:
	print("  Writing tree to", out)
for tip in T.get_terminals():
	name = re.sub("'", "", tip.name)
	tip.name = name
Phylo.write(T, out, "nexus")
subprocess.call(["sed", "-i", "-e", 's/\\\]//g', out])
subprocess.call(["sed", "-i", "-e", 's/\\[\\\//g', out])

if args.verbose:
	print("Done")
