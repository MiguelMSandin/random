#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Phylo
import os
import re

parser = argparse.ArgumentParser(description="Checks if all sequences from a fasta file are in the tree file and viceversa")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file in newick format.")

requiredArgs.add_argument("-f", "--fasta", dest="fasta", required=True,
						  help="A Fasta file.")

parser.add_argument("-c", "--check", dest="check", required=False, default=None, choices=['t', 'treeInFasta', 'f', 'fastaInTree'],
						  help="Only checks in one way: either if all tree tips are present in the fasta file ('t' or 'treeInFasta') or if all sequences in the fasta file are present in the tree file ('f' or 'fastaInTree').")

args = parser.parse_args()

# Reading files ------------------------------------------------------------------------------------
T = Phylo.read(args.tree, 'newick')
tree = list()
for line in T.get_terminals():
	tree.append(line.name)

fasta = list()
for line in SeqIO.parse(open(args.fasta), "fasta"):
	fasta.append(line.id)

# Start the search ---------------------------------------------------------------------------------
if args.check is None:
	checkT = True
	checkF = True
elif args.check == "t" or args.check == "treeInFasta":
	allGoodF2T = None
	checkT = True
	checkF = False
elif args.check == "f" or args.check == "fastaInTree":
	allGoodT2F = None
	checkT = False
	checkF = True

if checkT == True:
	allGoodT2F = True
	seqsT2F = set()
	for tip in tree:
		if tip not in fasta:
			allGoodT2F = False
			seqsT2F.add(tip)

if checkF == True:
	allGoodF2T = True
	seqsF2T = set()
	for seq in fasta:
		if seq not in tree:
			allGoodF2T = False
			seqsF2T.add(seq)

if allGoodT2F and allGoodF2T:
	print("All tree tips are found in the fasta file and viceversa.")
if allGoodT2F:
	print("All tree tips are found in the fasta file.")
if allGoodT2F == False:
	print("The following tree tips are not found in the fasta file:")
	if len(seqsT2F) < 10:
		for tip in seqsT2F:
			print("  ", tip)
	else:
		tmp = re.sub("\\.[^\\.]+$", "_tipsNotInFasta.list", args.tree)
		print("  More than 10 tips. Detailed file exported to:", tmp)
		with open(tmp, "w") as tmpo:
			for tip in seqsT2F:
				print(tip, file=tmpo)
if allGoodF2T:
	print("All sequences are found in the tree file.")
if allGoodF2T == False:
	print("The following sequences are not found in the tree file:")
	if len(seqsF2T) < 10:
		for seq in seqsF2T:
			print("  ", seq)
	else:
		tmp = re.sub("\\.[^\\.]+$", "_sequencesNotInTree.list", args.tree)
		print("  More than 10 sequences. Detailed file exported to:", tmp)
		with open(tmp, "w") as tmpo:
			for seq in seqsF2T:
				print(seq, file=tmpo)
print("Done")
