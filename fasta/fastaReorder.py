#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Phylo
import re
import sys

parser = argparse.ArgumentParser(description="Reorders a fasta file based on a tree or a list.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--fasta", dest="file_in", required=True,
					help="Fasta file.")

parser.add_argument("-t", "--tree", dest="tree", required=False, default=None,
					help="Tree file.")

parser.add_argument("-l", "--list", dest="list", required=False, default=None,
					help="A one column list.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
					help="The output fasta file name. By default will add '_treeOrder' or '_listOrder' to the fasta file name, depending on whether the ordering is performed after a tree or a list.")

parser.add_argument("-F", "--format", dest="format", required=False, default="newick",
					help="Tree file format (newick, nexus, nexml, phyloxml or cdao), default='newick'.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

if args.tree is None and args.list is None:
	print("\nError: Please select the ordering of the output file either from a phylogenetic tree ('-t/--tree') or a given list ('-l/--list').\nExiting\n")
	sys.exit(1)
if args.tree is not None and args.list is not None:
	print("\nError: Please select the ordering of the output file either from a phylogenetic tree ('-t/--tree') or a given list ('-l/--list'), but not both.\nExiting\n")
	sys.exit(1)

# Output file
if args.file_out is None:
	if args.tree is not None:
		outFile = re.sub("\\.[^\\.]+$", "_treeOrder.", args.file_in) + re.sub(".*\\.", "", args.file_in)
	if args.list is not None:
		outFile = re.sub("\\.[^\\.]+$", "_listOrder.", args.file_in) + re.sub(".*\\.", "", args.file_in)
else:
	outFile = args.file_out

# Reading the ordering file
ordering = list()
if args.tree is not None:
	if args.verbose:
		print("  Reading tree file")
	T = Phylo.read(args.tree, args.format)
	for line in T.get_terminals():
		tmp = line.name
		tmp = tmp.strip("'")
		ordering.append(line.name)
if args.list is not None:
	for line in open(args.list):
		tmp = line.strip("\n")
		ordering.append(tmp)

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
	for line in ordering:
		print(">" + str(line) + "\n" + str(fasta[line]), file=outfile)

if args.verbose:
	print("Done")
