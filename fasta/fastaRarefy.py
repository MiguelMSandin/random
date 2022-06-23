#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import statistics as st
import random

parser = argparse.ArgumentParser(description="From a fasta file, will export a table of rarefied observations. If an abundance table is given, abundances will be taken into consideration.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="fastaFile", required=True,
                    help="A fasta file.")

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

parser.add_argument("-i", "--identifier", dest="identifier", required=False, action="store_true",
                    help="If selected, will sample sequence identifiers and not sequences. Useful if (e.g.) there are duplicated sequences in your fasta file.")

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
	outFile = re.sub("\\.[^\\.]+$", "_rarefied.tsv", args.fastaFile)
else:
	outFile = args.output

rangeIter = [min(args.iterations), max(args.iterations)]
steps = list()
for i in range(rangeIter[0], rangeIter[1], args.steps):
	steps.append(i)
steps.append(rangeIter[1])

if args.identifier:
	reading = "identifiers"
else:
	reading = "sequences"

# Reading fasta ____________________________________________________________________________________
if args.verbose:
	print("  Reading fasta", args.fastaFile)
fasta = {}
for i in SeqIO.parse(open(args.fastaFile), "fasta"):
	fasta[i.id] = str(i.seq)

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
		print("  Replicating ", reading, " by abundance", end="")
		i = 0
		P = 0
	for key, value in fasta.items():
		if args.verbose:
			i += 1
			I = round(i/len(fasta)*100)
			if I > P:
				P = I
				print("\r  Replicating ", reading, " by abundance ", P, "%", sep="", end="")
		for j in range(int(abundance[key])):
			if args.identifier:
				reads.append(key)
			else:
				reads.append(value)
else:
	if args.verbose:
		print("  Extracting ", reading, end="")
		i = 0
		P = 0
	for key, value in fasta.items():
		if args.verbose:
			i += 1
			I = round(i/len(fasta)*100)
			if I > P:
				P = I
				print("\r  Extracting ", reading, " ", P, "%", sep="", end="")
		if args.identifier:
			reads.append(key)
		else:
			reads.append(value)
print("")

# Print information
if args.verbose:
	if args.abundance is not None:
		print("  Rarefication will be done:\n      -from", min(steps), "to", max(steps),
		"sampling size\n      -by steps of", steps[1]-steps[0],
		"\n      -with", args.replicates, "replicates\n      -in the total",
		len(reads), reading, "after replicating by abundance")
	else:
		print("  Rarefication will be done:\n      -from", min(steps), "to", max(steps),
		"sampling size\n      -by steps of", steps[1]-steps[0],
		"\n      -with", args.replicates, "replicates\n      -in the total",
		len(reads), reading)

# Test if it is possible to do not use replacement if selected
if args.replacement:
	if max(steps) > len(reads):
		print("\nWarning! You have selected a maximum sampling of", max(steps),
		"yet the sample has", len(reads),
		"reads.\nPlease consider using a smaller range or removing the replacement option.\nStopping\n")
		import sys
		sys.exit(1)

# Rarefying ________________________________________________________________________________________
out = {}
if args.verbose:
	print("  Rarefying", end="")
	i = 0
	P = 0
with open(outFile, 'w') as outfile:
	outfile.write("sampleSize\tmean\tsd\tmin\tp05\tp25\tp50\tp75\tp95\tmax\tcommon")
	for s in steps:
		if args.verbose:
			i += 1
			I = round(i/len(steps)*100)
			if I > P:
				P = I
				print("\r  Rarefying ", P, "%", sep="", end="")
		sample = list()
		common = set()
		for j in range(0, args.replicates):
			random.seed()
			if args.replacement:
				tmp = len(set(random.sample(reads, k=s)))
			else:
				tmp = len(set(random.choices(reads, k=s)))
			sample.append(tmp)
			common.add(tmp)
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
			 str(max(sample)) + '\t' +
			 str(len(common)))
		outfile.write(line)

if args.verbose:
	print("\n  Table exported to:", outFile)

# Summarising _____________________________________________________________________________________
if args.printSummary:
	print("  Final summary report:")
	print("    Fasta has a total of", len(fasta), "entries and", len(set(fasta.values())), "unique sequences")
	if args.abundance is not None:
		print("      (and", str(len(reads)), reading, "after replicating by abundance)")
	sampling = max([len(fasta), len(set(fasta.values()))])
	tmp = list()
	print("    Estimating", end="")
	for j in range(0, args.replicates):
		print("\r    Estimating ", j, "/", args.replicates, sep="", end="")
		if args.replacement:
			random.seed()
			tmp1 = len(set(random.sample(reads, k=sampling)))
		else:
			random.seed()
			tmp1 = len(set(random.choices(reads, k=sampling)))
		tmp.append(tmp1)
	tmp = st.mean(tmp)
	print("\r    When sampling ", sampling, " ", reading, ", an average of ", tmp, " unique ", reading, " are retrieved (", round(tmp/sampling*100, 2), "%)", sep="")
	if args.abundance is not None:
		tmp = list()
		print("    Estimating", end="")
		for j in range(0, args.replicates):
			print("\r    Estimating ", j, "/", args.replicates, sep="", end="")
			random.seed()
			if args.replacement:
				tmp1 = len(set(random.sample(reads, k=len(reads))))
			else:
				tmp1 = len(set(random.choices(reads, k=len(reads))))
			tmp.append(tmp1)
		tmp = st.mean(tmp)
		print("\r    When sampling ", len(reads), " ", reading, ", an average of ", tmp, " unique ", reading, " are retrieved (", round(tmp/sampling*100, 2), "%)", sep="")

# __________________________________________________________________________________________________
if args.verbose:
	print("Done")
