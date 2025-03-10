#!/usr/bin/env python3

import argparse
from Bio import Phylo

parser = argparse.ArgumentParser(description="Clusters a phylogenetic tree into subclades based on branch lengths, support, number of tips and/or relative branch lengths, ignoring the root node.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file in newick format.")

parser.add_argument("-b", "--branchLength", dest="branchLength", required=False, type=float, default=0,
					help="A minimum branch length. By default is 0.")

parser.add_argument("-c", "--confidence", dest="confidence", required=False, type=float, default=0,
					help="Minimum support values at nodes. By default is 0.")

parser.add_argument("-n", "--tips", dest="tips", required=False, default=1, type=int,
					help="Minimum number of tips in given clade. By default is 1.")

parser.add_argument("-r", "--relativeLength", dest="relativeLength", required=False, default=0, type=float,
					help="The relative branch length of a node regarding the average branch lengths within the node. By default is 0, which represents that a node has to have a branch length 0%% higher than the average branch lengths whithin the given node.")

# parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
# 					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

# parser.add_argument("-e", "--error", dest="error", required=False, action="store", type=int, default=None,
# 					help="If selected, will allow some error (in percentage: 0-100) for one threshold if the two other values pass thresholds.")

parser.add_argument("-o", "--output", dest="output", required=False, action="store", default=None,
					help="If selected, will export a list of nodes and their clades to the given output.")

parser.add_argument("-e", "--export", dest="export", required=False, action="store", default=None,
					help="If selected, will export tree with the selected nodes labelled to the given output.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Reading files ------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading tree")
# T = Phylo.read(args.tree, args.formaTree)
T = Phylo.read(args.tree, "newick")

# Checking if internal nodes pass the thresholds ---------------------------------------------------
subclades = 0
passLength = 0
passConfidence = 0
passTips = 0
passRelative = 0
pressenceLength = set()
pressenceConfidence = set()
internalNodes = 0
for clade in T.get_nonterminals():
	internalNodes += 1
	if internalNodes > 1:
		testLength = False
		testConfidence = False
		testTips = False
		testRelative = False
		avoidLength = False
		avoidConfidence = False
		avoidRelative = False
		if clade.branch_length is not None:
			pressenceLength.add(1)
			if clade.branch_length >= args.branchLength:
				testLength = True
				passLength += 1
		else:
			avoidLength = True
			pressenceLength.add(0)
		if clade.confidence is not None:
			pressenceConfidence.add(1)
			if clade.confidence >= args.confidence:
				testConfidence = True
				passConfidence += 1
		else:
			avoidConfidence = True
			pressenceConfidence.add(0)
		if clade.count_terminals() >= args.tips:
			testTips = True
			passTips += 1
		if args.relativeLength != 0 and clade.branch_length is not None:
			import statistics
			lengths = list()
			for sclade in clade.get_nonterminals():
				lengths.append(sclade.branch_length)
			for sclade in clade.get_terminals():
				lengths.append(sclade.branch_length)
			averageLength = statistics.mean(lengths)
			relAverageLength = averageLength + averageLength * args.relativeLength
			if clade.branch_length > relAverageLength:
				testRelative = True
				passRelative += 1
		else:
			avoidRelative = True
		if (testLength or avoidLength) and (testConfidence or avoidConfidence) and testTips and (testRelative or avoidRelative):
			subclades += 1
			if args.export is None:
				clade.name = "subclade" + str(subclades)
			else:
				clade.comment = '&!name="subclade_' + str(subclades) + '"'

# Exporting values ---------------------------------------------------------------------------------
print("  From", internalNodes, "internal nodes,", subclades, "passed all thresholds")
if sum(pressenceLength) == 1:
	print("    Nodes with branch lengths above ", args.branchLength, ": \t\t", passLength, sep="")
else:
	print("    Branch lengths not found in tree and therefore were ignored")
if sum(pressenceConfidence) == 1:
	print("    Nodes with confidence above ", args.confidence, ": \t\t\t", passConfidence, sep="")
else:
	print("    Confidence values not found in tree and therefore were ignored")
print("    Nodes with more than ", args.tips, " tips: \t\t\t", passTips, sep="")
if sum(pressenceLength) == 1 and args.relativeLength != 0:
	print("    Nodes with ", args.relativeLength*100, "% higher relative branch length: \t", passRelative, sep="")
if args.output is not None:
	with open(args.output, "w") as outfile:
		for clade in T.get_nonterminals():
			if args.export is None:
				if clade.name is not None:
					for tip in clade.get_terminals():
						print(str(clade.name) + "\t" + str(tip.name), file=outfile)
			else:
				if clade.comment is not None:
					for tip in clade.get_terminals():
						print(str(clade.name) + "\t" + str(tip.name), file=outfile)

# Writing file -------------------------------------------------------------------------------------
if args.export:
	if args.verbose:
		print("  Exporting labelled tree to", args.export)
	Phylo.write(T, args.export, "nexus")

if args.verbose:
	print("Done")
