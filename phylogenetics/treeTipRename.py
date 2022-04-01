#!/usr/bin/env python3

import argparse
from Bio import Phylo
import re

parser = argparse.ArgumentParser(description="")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

requiredArgs.add_argument("-l", "--list", dest="table", required=True,
					help="A tab delimited file with the tips to be renamed in one column and the new name in a second column.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-o", "--output", dest="output", required=False, default=None,
					help="The output file name. By default will add to the file name '_renamed' and the original extension.")

parser.add_argument("-F", "--formatOut", dest="formaTreeOut", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Setting variables --------------------------------------------------------------------------------

# Output file
if args.output is None:
	outFile = re.sub("\\.[^\\.]+$", "_renamed.", args.tree) + re.sub(".*\\.", ".", args.tree)
else:
	outFile = args.output

# Reading file -------------------------------------------------------------------------------------
if args.verbose:
	print("Reading tree file:", args.tree)
T = Phylo.read(args.tree, args.formaTree)

if args.verbose:
	print("Reading list file:", args.table)
names = {}
for line in open(args.table):
	tmp = line.strip().split()
	names[tmp[0]]=tmp[1]

# Rename tree --------------------------------------------------------------------------------------
if args.verbose:
	print("  Renaming")
notFound = list()
for line in T.get_terminals():
	tip = line.name
	tip = re.sub("\'", "", tip)
	if tip in names:
		line.name = names[tip]
	else:
		notFound.append(tip)

# Exporting renamed tree ---------------------------------------------------------------------------
if args.verbose:
	print("  Writing pruned tree to", outFile)
Phylo.write(T, outFile, args.formaTree)

# Checking errors and tips not foud ----------------------------------------------------------------
if len(notFound) > 0:
	print("  The following tips were not found in the list, and therefore not renamed:")
	for tip in notFound:
		print("  -", tip)

if args.verbose:
	print("Done")
