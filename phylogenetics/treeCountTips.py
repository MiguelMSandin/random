#!/usr/bin/env python3

import argparse
from Bio import Phylo

parser = argparse.ArgumentParser(description="Given a tree or a list of tree names in newick format, counts the number of tips and prints to the console.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", nargs='+', required=True,
						  help="A tree file(s).")

args = parser.parse_args()


for tree in args.tree:
	T = Phylo.read(tree, "newick")
	print(T.count_terminals(), "\t", tree)
