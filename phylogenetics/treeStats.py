#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Phylo
import statistics
import numpy as np

parser = argparse.ArgumentParser(description="Prints to the console some statistics of the tree")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--trees", dest="trees", required=True, nargs="+",
						  help="One or several tree files.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree(s) file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao. If more than one tree, all trees should have the same format")

parser.add_argument("-r", "--round", dest="round", required=False, default='16',
					help="If selected, will round the output values to the desired value. By default will show up to the 16th decimal point if available.")

parser.add_argument("-i", "--internal", dest="internal", required=False, action="store_false",
					help="If selected, will NOT take into account internal branches.")

parser.add_argument("-s", "--short", dest="short", required=False, action="store_true",
					help="If selected, will only print few comparable statistics.")

args = parser.parse_args()

if args.short:
	print("File\tTips\tNodes\tTotal branchlength")

if args.round == "16" and not args.short:
	r = 16
if args.round == "16" and args.short:
	r = 2
else:
	r = int(args.round)

# Reading file -------------------------------------------------------------------------------------
for tree in args.trees:
	T = Phylo.read(tree, args.formaTree)
	
	# Calculating stats --------------------------------------------------------------------------------
	# Number of tips
	ntips = T.count_terminals()
	
	# Number of nodes
	nnodes = len(T.get_nonterminals())
	
	# Branch lengths
	tips = list()
	lengths = list()
	for line in T.get_terminals():
		tips.append(line.name)
		lengths.append(line.branch_length)
	
	# Branch lengths of iternal branches
	if args.internal:
		for line in T.get_nonterminals():
			tmp = line.branch_length
			if tmp is not None:
				lengths.append(tmp)
		# Total branch length
		tbranchl = T.total_branch_length()
	else:
		tbranchl = 0
		for i in lengths:
			tbranchl = tbranchl + i
	
	# Printing information -----------------------------------------------------------------------------
	if args.short:
		print(tree, "\t", ntips, "\t", nnodes, "\t", round(tbranchl, r))
	else:
		print("File:\t", tree)
		
		count = 0
		if T.is_bifurcating():
			print("The tree is strictly bifurcating")
		else:
			print("The tree is NOT strictly bifurcating")
			for i in lengths:
				if i == 0:
					count += 1
		
		if T.rooted:
			print("The tree is rooted")
		else:
			print("The tree is NOT rooted")
		
		if args.internal == False:
			print("Internal branches not considered")
		
		print("")
		print("Number of tips:        \t", ntips)
		print("Number of nodes:       \t", nnodes)
		print("")
		print("Total branch lengths:  \t", round(tbranchl, r))
		print("")
		
		if None in lengths:
			print("Warning! Some (or all) branch lengths are not found")
		else:
			if count != 0:
				print("Shortest branch length:\t", round(min(lengths), r), "\t(found", count, "times)")
			else:
				print("Shortest branch length:\t", round(min(lengths), r))
			print("  5th percentile:      \t", round(np.percentile(lengths, 5), r))
			print("  25th percentile:     \t", round(np.percentile(lengths, 25), r))
			print("  Median:              \t", round(np.percentile(lengths, 50), r))
			print("Average branch length: \t", round(statistics.mean(lengths), r))
			print("  Standard deviation:  \t", round(statistics.stdev(lengths), r))
			print("  75th percentile:     \t", round(np.percentile(lengths, 75), r))
			print("  95th percentile:     \t", round(np.percentile(lengths, 95), r))
			print("Longest branch length: \t", round(max(lengths), r))
		print("")
