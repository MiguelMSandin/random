#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Phylo

parser = argparse.ArgumentParser(description="Prunes a tree from tips present in a list.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

requiredArgs.add_argument("-l", "--list", dest="list", required=True,
						  help="A list of tips to be removed from the tree (each tip should be given in a different line).")

parser.add_argument("-o", "--output", dest="out", required=False, action="store",
					help="The output tree file. If not selected, will add '_pruned.tre' before the extension of the input tree file.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick, nexus, nexml, phyloxml or cdao. Default='newick'")

parser.add_argument("-F", "--formatOut", dest="formatOut", required=False,
					help="The format of the output file: accepted formats are: newick, nexus, nexml, phyloxml or cdao. Default: will take input format.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_true",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# setting variables --------------------------------------------------------------------------------
if args.out is None:
	import re
	out = re.sub("\\.[^\\.]+$", "_pruned.tre", args.tree)
else:
	out = args.out

if args.formatOut is None:
	formatOut = args.formaTree
else:
	formatOut = args.formatOut
	
# Reading files ------------------------------------------------------------------------------------
if not args.verbose:
	print("  Reading files")

T = Phylo.read(args.tree, args.formaTree)
tips_in = T.count_terminals()

tips = [line.strip() for line in open(args.list)]
tips_list = len(tips)

# Prunning -----------------------------------------------------------------------------------------
if not args.verbose:
	print("  Prunning")
for tip in tips:
	T.prune(tip)

tips_out = T.count_terminals()

if not args.verbose:
	print("    Tips in: ", tips_in)
	print("    Tips out:", tips_out)
	print("    In total '", tips_list, "' tips were removed", sep="")
	print("  Writting file to: ", out)

# Writing file -------------------------------------------------------------------------------------
Phylo.write(T, out, formatOut)

if not args.verbose:
	print("Done")
