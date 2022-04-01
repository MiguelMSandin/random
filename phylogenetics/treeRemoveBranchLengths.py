#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Phylo

parser = argparse.ArgumentParser(description="Removes branch lengths of a phylogenetic tree.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

parser.add_argument("-o", "--output", dest="out", required=False, action="store",
					help="The output tree file. If not selected, will add '_noBrchL.tre' before the extension of the input tree file.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# setting variables --------------------------------------------------------------------------------
if args.out is None:
	import re
	out = re.sub("\\.[^\\.]+$", "_noBrchL", args.tree) + re.sub(".*\\.", ".", args.tree)
else:
	out = args.out

# Reading files ------------------------------------------------------------------------------------
if not args.verbose:
	print("  Reading files")

T = Phylo.read(args.tree, args.formaTree)

# Writing file -------------------------------------------------------------------------------------
if args.verbose:
	print("  Writting file to: ", out)

Phylo.write(T, out, args.formaTree, plain=True)

if args.verbose:
	print("Done")
