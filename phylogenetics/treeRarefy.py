#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Phylo
import statistics as st
import random

parser = argparse.ArgumentParser(description="From a tree file, will export a table of selected distances from rarefied tree tips. If an abundance table is given, abundances will be taken into consideration.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
                    help="A fasta file.")

parser.add_argument("-f", "--format", dest="formaTree", required=False,
					default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-a", "--abundance", dest="abundance", required=False, default=None,
                    help="If selected, will accomodate sequence abundance for the rarification. Then a tab separated table is needed with two columns: the name of the sequence and the abundance.")

parser.add_argument("-d", "--distance", dest="distance", required=False,
					default="root",
					choices=['root', 'branch', 'tips'],
                    help="The distance to be computed, either: 'root' (default) will calculate the total length from tip to root for every tip, 'branch' will calculate the branch length of every terminal tip, 'tips' will count the number of unique tips")

parser.add_argument("-o", "--output", dest="output", required=False, default=None,
                    help="The output name of the rarefied table. By default will add '_rarefiedDistance.tsv' to file name excluding the extension, being 'Distance' the chosen distance.")

parser.add_argument("-r", "--range", dest="iterations", required=False, type=int, nargs=2,
					default=[0, 20000],
                    help="The range of iterations. Will take the minimum and maximum values and create a range of 'steps' for random sampling. Default: '0 20000'.")

parser.add_argument("-s", "--steps", dest="steps", required=False, type=int,
					default="100",
                    help="The increase of sampling steps. Default: '100'.")

parser.add_argument("-n", "--replicates", dest="replicates", required=False, type=int,
					default="100",
                    help="The number of rarefaction replicates. Default: '100'.")

parser.add_argument("-m", "--normalize", dest="normalize", required=False, action="store_true",
                    help="If selected, the distance will be averaged by number of rarefied tips. This will have no effect for distance='tips'.")

parser.add_argument("-R", "--replacement", dest="replacement", required=False, action="store_true",
                    help="If selected, the random sampling will be done without replacement.")

parser.add_argument("-p", "--printSummary", dest="printSummary", required=False, action="store_false",
                    help="If selected, will not print a summary at the end.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
                    help="If selected, will not print information in the console.")

args = parser.parse_args()

# Setting parameters _______________________________________________________________________________
if args.output is None:
	if args.distance == "root":
		ext = "Root"
	elif args.distance == "branch":
		ext = "Branch"
	elif args.distance == "tips":
		ext = "Tips"
	else:
		ext = ""
	if args.normalize:
		n = "N"
	else:
		n = ""
	import re
	outFile = re.sub("\\.[^\\.]+$", "_rarefied", args.tree) + ext + n + ".tsv"
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

if args.distance != "tips":
	distances = {}
	if args.verbose:
		if args.distance == "root":
			print("  Getting branch legths to root", end="")
		if args.distance == "branch":
			print("  Getting branch legths", end="")
		i = 0
		P = 0
	for tip in T.get_terminals():
		if args.verbose:
				i += 1
				I = round(i/tipsN*100)
				if I > P:
					P = I
					if args.distance == "root":
						print("\r  Getting branch legths to root ", P, "%", sep="", end="")
					if args.distance == "branch":
						print("\r  Getting branch legths ", P, "%", sep="", end="")
		if args.distance == "root":
			distances[tip.name] = T.distance(tip)
		if args.distance == "branch":
			distances[tip.name] = tip.branch_length
	if args.verbose:
		print("")

if args.verbose:
	if args.distance == "root" or args.distance == "branch":
		length = list()
		for l in distances.values():
			length.append(l)
	if args.distance == "root":
		print("  The tree has", tipsN, "tips,", round(sum(length), 2), "total branch length from tip to root and", round(st.mean(length), 2), "on average.")
	if args.distance == "branch":
		print("  The tree has", tipsN, "tips,", round(sum(length), 2), "total branch lengths and", round(st.mean(length), 2), "on average.")
	if args.distance == "tips":
		print("  The tree has", tipsN, "tips")

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

# If selected replacement option, test if it is possible to use
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
	outfile.write("sampleSize\tmean\tsd\tmin\tp05\tp25\tp50\tp75\tp95\tmax\n")
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
			if args.distance == "tips":
				lengthOut = len(rarefied)
			else:
				for tip in rarefied:
					lengths.append(distances[tip])
				if args.normalize:
					if len(lengths) == 0:
						lengthOut = 0
					else:
						lengthOut = st.mean(lengths)
				else:
					lengthOut = sum(lengths)
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

# __________________________________________________________________________________________________
if args.verbose:
	print("Done")
