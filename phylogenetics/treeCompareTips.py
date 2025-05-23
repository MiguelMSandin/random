#!/usr/bin/env python3

import argparse
from Bio import Phylo

parser = argparse.ArgumentParser(description="Given a list of tree names in newick format, exports tips present in all tree or those not present in all trees.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", nargs='+', required=True,
						  help="A tree file(s).")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
                    help="Output list file. If not given will print tip names to console.")

parser.add_argument("-e", "--export", dest="export", required=False, default="rares",
					choices=['c', 'common', 'r', 'rares'],
                    help="Will export either all common tip names ('common' or 'c') to all trees or those not present in at least one tree ('rares' or 'r', default).")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Reading trees ____________________________________________________________________________________
if args.verbose:
	print("  Reading trees", end="")
tips = {}
i = 0
for tree in args.tree:
	i += 1
	if args.verbose:
		print("\r  Reading trees ", i, "/", len(args.tree), sep="", end="")
	T = Phylo.read(tree, "newick")
	for t in T.get_terminals():
		name = t.name
		if name not in tips.keys():
			tips[name] = 1
		else:
			tips[name] += 1

# Checking tips ____________________________________________________________________________________
if args.verbose:
	print("\n  Checking tips")
trees = len(args.tree)
out = set()
for key, value in tips.items():
	if args.export == "r" or args.export == "rares":
		if value < trees:
			out.add(key)
	if args.export == "c" or args.export == "common":
		if value == trees:
			out.add(key)

# Exporting ________________________________________________________________________________________
if args.file_out is None:
	if len(out) != 0:
		if args.verbose:
			if args.export == "r" or args.export == "rares":
				print("  The following tips do NOT appear in all trees:")
			if args.export == "c" or args.export == "common":
				print("  The following tips are common to all trees:")
		for tip in out:
			print(tip)
else:
	if args.verbose:
		if args.export == "r" or args.export == "rares":
			print("  Exporting tips that do NOT appear in all trees to", args.file_out)
		if args.export == "c" or args.export == "common":
			print("  Exporting tips that are common to all trees to", args.file_out)
	with open(args.file_out, "w") as outfile:
		for tip in out:
			print(tip, file=outfile)

if args.verbose:
	if len(out) == 0:
		print("  All tips are present in all trees", sep="")
	else:
		if args.export == "r" or args.export == "rares":
			print("  In total ", len(out), " tips do NOT appear in all trees", sep="")
		if args.export == "c" or args.export == "common":
			print("  In total ", len(out), " tips are common in all trees", sep="")
		print("    These represent ", round(len(out)/len(tips)*100, 2), "% of all unique tips", sep="")
	print("Done")
