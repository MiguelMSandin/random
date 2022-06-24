#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Phylo
import statistics as st
import random

parser = argparse.ArgumentParser(description="From a tree file, will export a table of branch lengths (from tip to root) from pruned trees after rarefying the tips. If an abundance table is given, abundances will be taken into consideration.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
                    help="A fasta file.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-a", "--abundance", dest="abundance", required=False, default=None,
                    help="If selected, will accomodate sequence abundance for the rarification. Then a tab separated table will be needed with two columns: the name of the sequence and the abundance.")

parser.add_argument("-o", "--output", dest="output", required=False, default=None,
                    help="The output name of the rarefied table. By default will add '_rarefied.tsv' to file name excluding the extension.")

parser.add_argument("-r", "--range", dest="iterations", required=False, type=int, nargs=2, default=[0, 20000],
                    help="The range of iterations. Will take the minimum and maximum values and create a range of 'steps' for random sampling. Default: '0 20000'.")

parser.add_argument("-s", "--steps", dest="steps", required=False, type=int, default="100",
                    help="The increase of sampling steps. Default: '100'.")

parser.add_argument("-n", "--replicates", dest="replicates", required=False, type=int, default="100",
                    help="The number of rarefaction replicates. Default: '100'.")

parser.add_argument("-N", "--normalize", dest="normalize", required=False, action="store_true",
                    help="If selected, the total branch length will be normalized by number of rarefied tips.")

parser.add_argument("-R", "--replacement", dest="replacement", required=False, action="store_true",
                    help="If selected, the random sampling will be done without replacement.")

parser.add_argument("-p", "--printSummary", dest="printSummary", required=False, action="store_false",
                    help="If selected, will not print a summary at the end.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
                    help="If selected, will not print information in the console.")

args = parser.parse_args()

# Setting parameters _______________________________________________________________________________
if args.output is None:
	import re
	outFile = re.sub("\\.[^\\.]+$", "_rarefied.tsv", args.tree)
else:
	outFile = args.output

rangeIter = [min(args.iterations), max(args.iterations)]
steps = list()
for i in range(rangeIter[0], rangeIter[1], args.steps):
	steps.append(i)
steps.append(rangeIter[1])

# Reading fasta ____________________________________________________________________________________
if args.verbose:
	print("  Reading tree file:", args.tree)
T = Phylo.read(args.tree, args.formaTree)
tipsN = T.count_terminals()

branchLengths = {}
if args.verbose:
	print("  Getting branch legths to root", end="")
	i = 0
	P = 0
for tip in T.get_terminals():
	if args.verbose:
			i += 1
			I = round(i/tipsN*100)
			if I > P:
				P = I
				print("\r  Getting branch legths to root ", P, "%", sep="", end="")
	branchLengths[tip.name] = T.distance(tip)
if args.verbose:
	print("")

# Reading abundances _______________________________________________________________________________
if args.abundance is not None:
	abundance = {}
	if args.verbose:
		print("  Reading abundance table", args.abundance)
	for line in open(args.abundance):
		line = line.strip().split()
		abundance[line[0]] = line[1]

# Extracting reads _________________________________________________________________________________
reads = list()
if args.abundance is not None:
	if args.verbose:
		print("  Replicating tips by abundance", end="")
		i = 0
		P = 0
	for tip in T.get_terminals():
		if args.verbose:
			i += 1
			I = round(i/tipsN*100)
			if I > P:
				P = I
				print("\r  Replicating tips by abundance ", P, "%", sep="", end="")
		for j in range(int(abundance[tip.name])):
			reads.append(tip.name)
else:
	if args.verbose:
		print("  Extracting tips", end="")
		i = 0
		P = 0
	for tip in T.get_terminals():
		if args.verbose:
			i += 1
			I = round(i/tipsN*100)
			if I > P:
				P = I
				print("\r  Extracting tips ", P, "%", sep="", end="")
			reads.append(tip.name)
if args.verbose:
	print("")

# Print information
if args.verbose:
	if args.abundance is not None:
		print("  Rarefication will be done:\n      -from", min(steps), "to", max(steps),
		"sampling size\n      -by steps of", steps[1]-steps[0],
		"\n      -with", args.replicates, "replicates\n      -in the total",
		len(reads), "tips after replicating by abundance")
	else:
		print("  Rarefication will be done:\n      -from", min(steps), "to", max(steps),
		"sampling size\n      -by steps of", steps[1]-steps[0],
		"\n      -with", args.replicates, "replicates\n      -in the total",
		len(reads), "tips")

# Test if it is possible to do not use replacement if selected
if args.replacement:
	if max(steps) > len(reads):
		print("\nWarning! You have selected a maximum sampling of", max(steps),
		"yet the sample has", len(reads),
		"reads.\nPlease consider using a smaller range or removing the replacement option.\nStopping\n")
		import sys
		sys.exit(1)

# Rarefying ________________________________________________________________________________________
if args.verbose:
	print("  Rarefying", end="")
	i = 0
	P = 0
with open(outFile, 'w') as outfile:
	if args.normalize:
		outfile.write("sampleSize\tmeanNorm\tsdNorm\tmin\tp05\tp25\tp50\tp75\tp95\tmax\n")
	else:
		outfile.write("sampleSize\tmeanTotal\tsdTotal\tmin\tp05\tp25\tp50\tp75\tp95\tmax\n")
	for s in steps:
		if args.verbose:
			i += 1
			I = round(i/len(steps)*100)
			if I > P:
				P = I
				print("\r  Rarefying ", P, "%", sep="", end="")
		sample = list()
		for j in range(0, args.replicates):
			lengths = list()
			random.seed()
			if args.replacement:
				rarefied = set(random.sample(reads, k=s))
			else:
				rarefied = set(random.choices(reads, k=s))
			for tip in rarefied:
				lengths.append(branchLengths[tip])
			lengthOut = sum(lengths)
			if args.normalize:
				lengthOut = lengthOut/s
			sample.append(lengthOut)
		sort = sorted(sample)
		line = str(str(s) + '\t' +
			 str(st.mean(sample)) + '\t' +
			 str(st.stdev(sample)) + '\t' +
			 str(min(sample)) + '\t' +
			 str(sort[int(len(sample)*0.05)]) + '\t' +
			 str(sort[int(len(sample)*0.25)]) + '\t' +
			 str(sort[int(len(sample)*0.5)]) + '\t' +
			 str(sort[int(len(sample)*0.75)]) + '\t' +
			 str(sort[int(len(sample)*0.95)]) + '\t' +
			 str(max(sample)) + '\n')
		outfile.write(line)

if args.verbose:
	print("\n  Table exported to:", outFile)

if args.printSummary:
	print("  Final summary report:")
	total = 0
	branchLengths
	for l in branchLengths.values():
		total += l
	print("    Tree has a total of", tipsN, "tips and", round(total, 2), "total branch length")
	sample = list()
	for j in range(0, args.replicates):
		random.seed()
		if args.replacement:
			rarefied = set(random.sample(reads, k=tipsN))
		else:
			rarefied = set(random.choices(reads, k=tipsN))
		lengths = list()
		for tip in rarefied:
			lengths.append(branchLengths[tip])
		lengthOut = sum(lengths)
		if args.normalize:
			lengthOut = lengthOut/s
		sample.append(lengthOut)
	tmp = st.mean(sample)
	print("\r    When sampling ", tipsN, " tips, an average branch length of ", round(tmp, 2), " is retrieved (", round(tmp/total*100, 2), "%)", sep="")
	if args.abundance is not None:
		sample = list()
		for j in range(0, args.replicates):
			random.seed()
			if args.replacement:
				rarefied = set(random.sample(reads, k=len(reads)))
			else:
				rarefied = set(random.choices(reads, k=len(reads)))
			lengths = list()
			for tip in rarefied:
				lengths.append(branchLengths[tip])
			lengthOut = sum(lengths)
			if args.normalize:
				lengthOut = lengthOut/s
			sample.append(lengthOut)
		tmp = st.mean(sample)
		print("\r    When sampling ", len(reads), " tips, an average branch length of ", round(tmp, 2), " is retrieved (", round(tmp/total*100, 2), "%)", sep="")

# __________________________________________________________________________________________________
if args.verbose:
	print("Done")
