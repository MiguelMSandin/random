#!/usr/bin/env python3

import argparse
from Bio import Phylo
import statistics
import math
import re

parser = argparse.ArgumentParser(description="From a time calibrated tree, it will export a tab delimited table with the number of Lineages Through Time (LTT).")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file in newick format.")

parser.add_argument("-o", "--out", dest="out", required=False, default=None,
						  help="The output file name. By default, will remove the extension of the input tree file and add '_LTT.tsv'.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default) and nexus.")

parser.add_argument("-s", "--subtrees", dest="subtrees", required=False, action="store_true",
					help="If selected, will also append to the table the LTT data for every annotated clade.")

parser.add_argument("-d", "--hpd", dest="hpd", required=False, action="store_true",
					help="If selected, will also append to the table the Highest Posterior Density probability for every time point.")

parser.add_argument("-u", "--ultrametric", dest="ultrametric", required=False, action="store_false",
					help="If selected, will not take the average of all distances from clade to tip but return an error if these are different.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# setting variables --------------------------------------------------------------------------------
if args.out is None:
	out = re.sub("\\.[^\\.]+$", "_LTT.tsv", args.tree)
else:
	out = args.out

# Reading files ------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading files")
T = Phylo.read(args.tree, args.formaTree)

# Extracting node ages -----------------------------------------------------------------------------
if args.verbose:
	print("  Extracting node ages")
if args.hpd:
	heights = {}
else:
	heights = list()
for clade in T.get_nonterminals():
	height = set()
	for tip in clade.get_terminals():
		dist = clade.distance(clade, tip)
		dist = round(dist, 6)
		height.add(dist)
	if len(height) != 1 and not args.ultrametric:
		import sys
		print("  Error: Different distances from node to tips. Please check your tree is ultrametric.")
		sys.exit(1)
	elif args.ultrametric:
		height = statistics.mean(height)
	else:
		height = next(iter(height))
	height = round(height, 6)
	if args.hpd:
		tmp = clade.comment
		tmp = re.sub(".*HPD=", "", tmp)
		tmp = re.sub("{", "", tmp)
		tmp = re.sub("}.*", "", tmp)
		hpd05 = re.sub(",.*", "", tmp)
		hpd95 = re.sub(".*,", "", tmp)
		heights[height] = [hpd05, hpd95]
	else:
		heights.append(height)

# Extracting subtrees ------------------------------------------------------------------------------
if args.subtrees:
	if args.verbose:
		print("  Extracting subtrees")
	subtrees = {}
	for clade in T.get_nonterminals():
		if clade.comment is not None:
			comment = clade.comment
			if "name" in comment:
				name = re.sub('.*name="', "", comment)
				name = re.sub(',.*', "", name)
				name = re.sub('"', "", name)
				name = re.sub('\\]', "", name)
				if args.hpd:
					subtrees[name] = {}
				else:
					subtrees[name] = list()
				for subclade in clade.get_nonterminals():
					height = set()
					for tip in subclade.get_terminals():
						dist = subclade.distance(subclade, tip)
						dist = round(dist, 6)
						height.add(dist)
					height = statistics.mean(height)
					height = round(height, 6)
					if args.hpd:
						tmp = subclade.comment
						tmp = re.sub(".*HPD=", "", tmp)
						tmp = re.sub("{", "", tmp)
						tmp = re.sub("}.*", "", tmp)
						hpd05 = re.sub(",.*", "", tmp)
						hpd95 = re.sub(".*,", "", tmp)
						subtrees[name][height] = [hpd05, hpd95]
					else:
						subtrees[name].append(height)

# Writing file -------------------------------------------------------------------------------------
if args.verbose:
	print("  Writing table to", out)
e = math.exp(1)
with open(out, "w") as outfile:
	if not args.subtrees and args.hpd:
		print("time\tlineages\tlnLineages\thpd05\thpd95", file=outfile)
		n = 0
		for h in sorted(list(heights.keys()), reverse=True):
			n += 1
			ln = math.log(n, e)
			print(str(h) + "\t" + str(n) + "\t" + str(ln) + "\t" + str(heights[h][0]) + "\t" + str(heights[h][1]), file=outfile)
		n += 1
		ln = math.log(n, e)
		print("0\t" + str(n) + "\t" + str(ln) + "\t0\t0", file=outfile)
	if not args.subtrees and not args.hpd:
		print("time\tlineages\tlnLineages", file=outfile)
		n = 0
		for h in sorted(heights, reverse=True):
			n += 1
			ln = math.log(n, e)
			print(str(h) + "\t" + str(n) + "\t" + str(ln), file=outfile)
		n += 1
		ln = math.log(n, e)
		print("0\t" + str(n) + "\t" + str(ln), file=outfile)
	if args.subtrees and args.hpd:
		print("tree\ttime\tlineages\tlnLineages\thpd05\thpd95", file=outfile)
		n = 0
		for h in sorted(list(heights.keys()), reverse=True):
			n += 1
			ln = math.log(n, e)
			print("main\t" + str(h) + "\t" + str(n) + "\t" + str(ln) + "\t" + str(heights[h][0]) + "\t" + str(heights[h][1]), file=outfile)
		n += 1
		ln = math.log(n, e)
		print("main\t0\t" + str(n) + "\t" + str(ln) + "\t0\t0", file=outfile)
		for key, dicts in subtrees.items():
			n = 0
			for h in sorted(list(dicts.keys()), reverse=True):
				n += 1
				ln = math.log(n, e)
				value = subtrees[key][h]
				print(str(key) + "\t" + str(h) + "\t" + str(n) + "\t" + str(ln) + "\t" + str(value[0]) + "\t" + str(value[1]), file=outfile)
			n += 1
			ln = math.log(n, e)
			print(str(key) + "\t0\t" + str(n) + "\t" + str(ln) + "\t0\t0", file=outfile)
	if args.subtrees and not args.hpd:
		print("tree\ttime\tlineages\tlnLineages", file=outfile)
		n = 0
		for h in sorted(heights, reverse=True):
			n += 1
			ln = math.log(n, e)
			print("main\t" + str(h) + "\t" + str(n) + "\t" + str(ln), file=outfile)
		for key, values in subtrees.items():
			n = 0
			for h in sorted(values, reverse=True):
				n += 1
				ln = math.log(n, e)
				print(str(key) + "\t" + str(h) + "\t" + str(n) + "\t" + str(ln), file=outfile)
			n += 1
			ln = math.log(n, e)
			print(str(key) + "\t0\t" + str(n) + "\t" + str(ln), file=outfile)
if args.verbose:
	print("Done")
