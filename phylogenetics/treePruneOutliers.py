#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Phylo
import statistics
import numpy as np
import sys
import re

parser = argparse.ArgumentParser(description="Prunes a tree of tips which their branch length are identified as outliers by either the Z-scores, the interquartile range or the generalized extreme studentized deviate method.",
								 epilog="*Depending on the method, outliers are defined if; Z-scores:|i—μ|/σ > t; IQR: i < q1-(t*iqr) OR i > q3+(t*iqr) (being 'i' the given branch length, 'µ' the average branch length, 'σ' the standard deviation, 't' the chosen threshold, 'q1' the 25th quartile, 'q3' the 75th quartile and 'iqr' the difference between 'q3' and 'q1'); and gESD: Rosner, Bernard (1983), Percentage Points for a Generalized ESD Many-Outlier Procedure,Technometrics, 25(2), pp. 165-172.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-o", "--output", dest="output", required=False, default=None,
					help="The output file name. By default will add '_pruned' to the file name. If output='false', the tree will not be exported.")

parser.add_argument("-m", "--method", dest="method", required=False, default='zscores', choices=['zscores', 'iqr', 'gesd'],
					help="The method to identify outliers: Either by Z-scores ('zscores', default), by the InterQuartile Range ('iqr') or by the generalized Extreme Studentized Deviate (gESD: 'gesd')")

parser.add_argument("-r", "--threshold", dest="threshold", required=False, default=None,
					help="The threshold to identify outliers. By default: zscores=2; iqr:1.5; gesd: 0.05 (refers to the significance)")

parser.add_argument("-n", "--maximum", dest="maximum", required=False, default=None,
					help="For gESD method only, an estimate of the maximum number of outliers in the dataset. Default= 1/10 of the number of tips.")

parser.add_argument("-l", "--list", dest="table", required=False, action="store_true",
					help="If selected, will write the tips, its branch lengths and whether it was considered an outlier or not a tab delimited table to the input file adding '_branchLengths.tsv'")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

methodsHelp = parser.add_argument_group('Methods')

args = parser.parse_args()

# Reading file -------------------------------------------------------------------------------------
if args.verbose:
	print("Reading tree file:", args.tree)
T = Phylo.read(args.tree, args.formaTree)

# Extracting branch lengths ------------------------------------------------------------------------
if args.verbose:
	print("  Getting branch lengths")
tips = list()
lengths = list()
for line in T.get_terminals():
	tips.append(line.name)
	lengths.append(line.branch_length)

# Setting variables --------------------------------------------------------------------------------

# Output file
if args.output is None:
        outFile = re.sub("\\.[^\\.]+$", "_pruned.", args.tree) + re.sub(".*\\.", "", args.tree)
else:
        outFile = args.output

# Branch lengths file if applicable
if args.table:
        outTable = re.sub("\\.[^\\.]+$", "_branchLengths.tsv", args.tree)
else:
        outTable = args.table

# Z-Scores
if args.method == "zscores":
	a = statistics.mean(lengths)
	s = statistics.stdev(lengths)
	if args.threshold is None:
		t = 2
		ext = "(default value)"
	else:
		t = float(args.threshold)
		ext = ""
	if args.verbose:
		print("  Using Z-scores method for outlier identification")
		print("    Average:           ", a)
		print("    Standard deviation:", s)
		print("    Threshold:         ", t, ext)
		print("      Formula: | i —", round(a, 2), " | /", round(s, 2), " >", round(t, 2))

# IQR
if args.method == "iqr":
	
	q1 = np.percentile(lengths, 25)
	q3 = np.percentile(lengths, 75)
	iqr = q3 - q1
	if args.threshold is None:
		t = 1.5
		ext = "(default value)"
	else:
		t = float(args.threshold)
		ext =""
	lower = q1-(t*iqr)
	upper = q3+(t*iqr)
	if args.verbose:
		print("  Using IQR method for outlier identification")
		print("    25th quartile:", q1)
		print("    75th quartile:", q3)
		print("    Threshold:    ", t, ext)
		print("      Formula: i <", round(q1, 4), "- (", t, "*", round(iqr, 4), ") OR i >", round(q3, 4),"+ (", t, "*", round(iqr, 4), ") = i <", round(lower, 4), "| i >", round(upper, 4))

# gESD
if args.method == "gesd":
	from PyAstronomy import pyasl
	if args.threshold is None:
		t = 0.05
		ext1 = "(default value)"
	else:
		t = float(args.threshold)
		ext1 = ""
	if args.maximum is None:
		m = int(len(tips)/10)
		ext2 = "(default value: 10%)"
	else:
		m = int(args.maximum)
		ext2 = ""
	if args.verbose:
		print("  Using gESD method for outlier identification")
		print("    Significance:                      ", t, ext1)
		print("    Maximum number of outliers allowed:", m, ext2)


# Calculating outliers -----------------------------------------------------------------------------
toPrune = list()
zips = zip(tips, lengths)

# Z-scores
if args.method == "zscores":
	for tip, length in zips:
		j = abs((length - a)) / s
		if j > t:
			toPrune.append(tip)

# IQR
if args.method == "iqr":
	for tip, length in zips:
		if length < lower or length > upper:
			toPrune.append(tip)

# gESD
if args.method == "gesd":
	gesd = pyasl.generalizedESD(lengths, m, t)
	for j in gesd[1]:
		toPrune.append(tips[j])

# Pruning ------------------------------------------------------------------------------------------
if outFile == "false":
	if args.verbose:
		print("  In total", len(toPrune), "tips were considered as outliers")
if outFile != "false":
	if args.verbose:
		print("  Prunning a total of", len(toPrune), "tips...")
	for tip in toPrune:
		T.prune(tip)

# Writing files ------------------------------------------------------------------------------------
if outFile == "false":
	if args.verbose:
		print("  Pruned tree is NOT exported")
if outFile != "false":
	if args.verbose:
		print("  Writing pruned tree to", outFile)
	Phylo.write(T, outFile, args.formaTree)

if args.table:
	outTable = re.sub("\\.[^\\.]+$", "_branchLengths.tsv", args.tree)
	if args.verbose:
		print("  Writing branch lengths to", outTable)
	identified = list()
	for tip in tips:
		if tip in toPrune:
			identified.append("Outlier")
		else:
			identified.append("")
	with open(outTable, "w") as outtable:
		zips = zip(tips, lengths, identified)
		for tip, length, iden in zips:
			print(str(tip) + '\t' + str(length) + '\t' + str(iden), file=outtable)

if args.verbose:
	if outFile != "false":
		print("  Tips in input tree: ", len(tips))
		print("  Tips in output tree:", T.count_terminals())
	print("Done")
