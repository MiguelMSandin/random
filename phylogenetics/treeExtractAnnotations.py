#!/usr/bin/env python3

import argparse
from Bio import Phylo
import re

parser = argparse.ArgumentParser(description="Extract a tsv table with the tip names in one column and all annotations until the root (separated by '|') of every given tip in a second column.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='nexus',
					help="The annotated tree file format: accepted formats are: nexus (default), newick, nexml, phyloxml or cdao.")

parser.add_argument("-o", "--output", dest="output", required=False, default=None,
					help="The output name of the annotated tree file. By default will take the name of the tree to annotate followed by '_annotateions.tsv'")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Reading files ------------------------------------------------------------------------------------
if args.verbose:
	print("Reading tree:  ", args.tree)
T = Phylo.read(args.tree, args.formaTree)

# Output file
if args.output is None:
	outFile = re.sub("\\.[^\\.]+$", "_annotations.tsv", args.tree)
else:
	outFile = args.output

# Loop through internal nodes and look for childs -------------------------------------------------
if args.verbose:
	print("  Extracting annotations and tip names to:", outFile)
with open(outFile, 'w') as outfile:
	for tip in T.get_terminals():
		label = tip.name
		nodes = T.get_path(tip)
		annotation = ""
		for node in nodes[:-1]:
			if "name" in node.comment:
				tmp = node.comment
				tmp = re.sub('.*name="', "", tmp)
				tmp = re.sub('".*', "", tmp)
			else:
				tmp= ""
			if annotation == "" and tmp != "":
				annotation = tmp
			elif annotation != "" and tmp != "":
				annotation = annotation + '|' + tmp
		print(label + "\t" + annotation, file=outfile)

if args.verbose:
	print("Done")
