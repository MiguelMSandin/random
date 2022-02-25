#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Phylo
import os
import re

parser = argparse.ArgumentParser(description="From every tree tip in a newick tree, searches every sequence name in a fasta file and exports a polytomic tree with the sequence names.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

requiredArgs.add_argument("-f", "--fasta", dest="fasta", required=True,
						  help="A Fasta file.")

parser.add_argument("-d", "--outputDir", dest="dir_out", required=False, action="store",
					help="The output direcotry file where it will export a text file for every tree tip containing the sequence names found in the fasta file.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# setting variables --------------------------------------------------------------------------------
if args.verbose:
	print("  Setting variables")

if args.dir_out is None:
	import re
	outDir = re.sub("\\.[^\\.]+$", "", args.tree)
else:
	outDir = args.dir_out

if os.path.exists(outDir):
	print("    Output directory already exists, deleting...")
	import shutil
	shutil.rmtree(outDir)
	os.mkdir(outDir)
else:
	os.mkdir(outDir)

# Reading files ------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading tree")
T = Phylo.read(args.tree, 'newick')

if args.verbose:
	print("  Reading fasta")
fasta = list()
for line in SeqIO.parse(open(args.fasta), "fasta"):
	fasta.append(line.id)

# Start the search ---------------------------------------------------------------------------------
if args.verbose:
	print("  Extracting names")
allGood = True
seqs = set()
excluded = list()
for line in T.get_terminals():
	sp = line.name
	warnings = 0
	repeated = 0
	seqsRep = set()
	tmpFile = outDir + "/" + sp + ".txt"
	with open(tmpFile, "w") as tmp:
		for i in fasta:
			match = re.search(sp, i)
			if match:
				if i in seqs:
					repeated += 1
					seqsRep.add(i)
				seqs.add(i)
				print(i, file=tmp)
			if match is None:
				warnings += 1
		if warnings == len(fasta):
			print("    Warning! Taxa '", sp, "' has no matching sequences", sep="")
			excluded.append(sp)
			allGood = False
		if repeated > 0:
			print("    Warning!! There are '", repeated, "' repeated sequences when looking for '", sp, "', also matching:", sep="")
			allGood = False
			names = set()
			for i in T.get_terminals():
				name = i.name
				string = re.compile(name)
				if name != sp:
					for j in seqsRep:
						if string.search(j):
							names.add(name)
			excluded.append(str(str(sp)+" also matches "+str(names)))
			print("     ", *names)

if allGood:
	print("")
	print("All tree tips are found in the fasta file and there are no repeated tree tips.")
	print("Ready to replace names :) !")
	print("")
else:
	excluded_out = outDir + "_notFound.list"
	with open(excluded_out, "w") as tmp:
		for sp in excluded:
			print(sp, file=tmp)
	print("")
	print("  Taxa not found in fasta file or conflicting with other taxa are exported to:", excluded_out)
	print("")

if args.verbose:
	print("Done")
