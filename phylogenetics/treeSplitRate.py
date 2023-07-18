#!/usr/bin/env python3

import argparse
from Bio import Phylo
import re

parser = argparse.ArgumentParser(description="From a time calibrated tree, it will export a tab delimited table with the split rate at every given time interval.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

parser.add_argument("-o", "--out", dest="out", required=False, default=None,
						  help="The output file name. By default, will remove the extension of the input tree file and add '_splitRateINT.tsv', where INT is the chosen time interval.")

parser.add_argument("-s", "--steps", dest="steps", required=False, type=int, default=50,
					help="The number of steps to estimate the split rate. Default=50 steps.")

parser.add_argument("-i", "--interval", dest="interval", required=False, type=int, default=None,
					help="The time interval within the split rate will be estimated. If used, will ignore the inpit given in '-s/--steps'.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default) and nexus.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Reading files ------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading file")
T = Phylo.read(args.tree, args.formaTree)

# Extracting node ages -----------------------------------------------------------------------------
if args.verbose:
	print("  Extracting node ages")
heights = list()
for clade in T.get_nonterminals():
	height = None
	for tip in clade.get_terminals():
		if height is None:
			height = clade.distance(clade, tip)
	heights.append(height)

i = T.count_terminals()
lineages = {}
for h in sorted(heights):
	lineages[h] = i
	i = i - 1

# setting variables --------------------------------------------------------------------------------
rootAge = max(heights)
if args.interval is None:
	inter = int(rootAge / args.steps)
	print("  Estimate split rate for", args.steps, "steps every", inter, "units of time")
else:
	inter = int(args.interval)
	print("  Estimate split rate every", inter, "units of time for a total of", round(rootAge/inter), "steps")

if args.out is None:
	out = re.sub("\\.[^\\.]+$", "_splitRate", args.tree)
	out = out + str(inter) + ".tsv"
else:
	out = args.out

# Writing file -------------------------------------------------------------------------------------
if args.verbose:
	print("  Writing table to", out)
with open(out, "w") as outfile:
	print("time\tlineages\tsplits\trate", file=outfile)
	for i in range(round(min(heights)), round(max(heights)), inter):
		if i > 0:
			j = i - inter
			l = 0
			s = 0
			for height, lineage in lineages.items():
				if height < i and height > j:
					s += 1
					if lineage > l:
						l = lineage
			if l == 0:
				for h2, l2 in lineages.items():
					if h2 > j:
						if l2 > l:
							l = l2
			print(str(i), "\t", str(l), "\t", str(s), "\t", str(s/l), file=outfile)

if args.verbose:
	print("Done")
