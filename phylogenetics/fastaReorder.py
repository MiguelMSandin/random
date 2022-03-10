#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Phylo
import re

parser = argparse.ArgumentParser(description="Reorders a fasta file based on a tree.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--fasta", dest="file_in", required=True,
					help="Fasta file.")

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
					help="Tree file.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
					help="The output fasta file name. By default will add '_treeOrder' to the fasta file name.")

parser.add_argument("-F", "--format", dest="format", required=False, default="newick",
					help="Tree file format (newick, nexus, nexml, phyloxml or cdao), default='newick'.")

parser.add_argument("-d", "--drawTree", dest="draw", required=False, action="store_true", default=None,
					help="If selected, draws the tree in the console while reordering the fasta file.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Output file
if args.file_out is None:
	outFile = re.sub("\\.[^\\.]+$", "_treeOrder.", args.file_in) + re.sub(".*\\.", "", args.file_in)
else:
	outFile = args.file_out

# Reading tree file
if args.verbose:
	print("  Reading tree file")
T = Phylo.read(args.tree, args.format)
if args.draw is not None:
	Phylo.draw_ascii(T)

# Reading fasta file
if args.verbose:
	print("  Reading fasta file")
fasta = {}
for line in SeqIO.parse(open(args.file_in), "fasta"):
	fasta[line.id] = str(line.seq)

# Reordering and exporting
if args.verbose:
	print("  Reordering and writing fasta file to:", outFile)
with open(outFile, "w") as outfile:
	for line in  T.get_terminals():
		sp = line.name
		sp = sp.strip("'")
		print(">" + str(sp) + "\n" + str(fasta[sp]), file=outfile)

if args.verbose:
	print("Done")
