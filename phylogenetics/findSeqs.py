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
allGoodT2F = True
seqsT2F = set()
for tip in tree:
	if tip not in fasta:
		allGoodT2F = False
		seqsT2F.add(tip)

allGoodF2T = True
seqsF2T = set()
for seq in fasta:
	if seq not in tree:
		allGoodF2T = False
		seqsF2T.add(seq)

if allGoodT2F and allGoodF2T:
	print("")
	print("All tree tips are found in the fasta file and viceversa.")
	print("")
if allGoodT2F == False:
	print("")
	print("The following tree tips are not found in the fasta file:")
	if len(seqsT2F) < 10:
		print("  ", *seqsT2F)
	else:
		tmp = re.sub("\\.[^\\.]+$", "_tipsNotInFasta.list", args.tree)
		print("  More than 10 tips. Detailed file exported to:", tmp)
		with open(tmp, "w") as tmpo:
			for tip in seqsT2F:
				print(tip, file=tmpo)
			
	print("")
if allGoodF2T == False:
	print("")
	print("The following sequences are not found in the tree file:")
	print("  ", *seqsF2T)
	if len(seqsF2T) < 10:
		print("  ", *seqsF2T)
	else:
		tmp = re.sub("\\.[^\\.]+$", "_sequencesNotInTree.list", args.tree)
		print("  More than 10 sequences. Detailed file exported to:", tmp)
		with open(tmp, "w") as tmpo:
			for seq in seqsF2T:
				print(seq, file=tmpo)
	print("")

