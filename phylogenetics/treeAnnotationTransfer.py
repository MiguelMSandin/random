#!/usr/bin/env python3

import argparse
from Bio import Phylo
import re

parser = argparse.ArgumentParser(description="Given an annotated tree, will take the annotations and transfer them to a new tree. More precisely, for every annotated node, will take the first an last tip name, look for the last common ancestor of these two tips in the tree to be annotated and copy the annotation. Therefore the two trees must have identical tip names (although not necessarily all the tips).")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="treeRef", required=True,
						  help="A tree file.")

requiredArgs.add_argument("-a", "--annotate", dest="treeToAnnot", required=True,
						  help="A tree file to be annotated.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The annotated tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-F", "--formatToAnnotate", dest="formaTreeToAnnotate", required=False, default='newick',
					help="The tree to annotate file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-o", "--output", dest="output", required=False, default=None,
					help="The output name of the annotated tree file. By default will take the name of the tree to annotate followed by '_annotated' and the extension.")

parser.add_argument("-O", "--formatOutput", dest="formatOutput", required=False, default='newick',
					help="The tree output format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

methodsHelp = parser.add_argument_group('Methods')

args = parser.parse_args()

# Reading files ------------------------------------------------------------------------------------
if args.verbose:
	print("Reading annotated tree:  ", args.treeRef)
T = Phylo.read(args.treeRef, args.formaTree)

if args.verbose:
	print("Reading tree to annotate:", args.treeToAnnot)
A = Phylo.read(args.treeToAnnot, args.formaTreeToAnnotate)

# Setting variables --------------------------------------------------------------------------------

# Output file
if args.output is None:
	outFile = re.sub("\\.[^\\.]+$", "_annotated.", args.treeToAnnot) + re.sub(".*\\.", "", args.treeToAnnot)
else:
	outFile = args.output

# Loop through internal nodes and look for childs -------------------------------------------------
if args.verbose:
	print("  Extracting annotations and tip names")

annot={}
for line in T.get_nonterminals():
	if line.name is not None:
		tips = line.count_terminals()
		tmp = list()
		for tip in line.get_terminals():
			tmp.append(tip.name)
		annot[line.name] = list(tmp)
	elif line.comment is not None:
		annotation = line.comment
		annotation = re.sub('\[|\]|"', "", annotation)
		annotation = annotation.split(",")[0]
		annotation = re.sub(".*=", "", annotation)
		if 'collapsed' not in annotation and 'rotate' not in line.comment:
			tips = line.count_terminals()
			tmp = list()
			for tip in line.get_terminals():
				tmp.append(re.sub("\'", "", tip.name))
			annot[annotation] = list(tmp)

leafs = list()
for line in A.get_terminals():
	leafs.append(line.name)

# Transfer the annotations ------------------------------------------------------------------------
if args.verbose:
	print("  Transferring annotations")

for annotation, tips in annot.items():
	first = None
	for tip in tips:
		if tip in leafs:
			first = tip
			break
	last = None
	for tip in reversed(tips):
		if tip in leafs:
			last = tip
			break
	try:
		clade = A.common_ancestor(first, last)
		clade.name = annotation
	except:
		pass

# Counting the annotations transferred ------------------------------------------------------------
annotated = list()
for line in A.get_nonterminals():
	name = line.name
	if name is not None:
		annotated.append(name)

# Finding not transferred annotations --------------------------------------------------------------
notFound = list()
if len(annotated) != len(annot):
	for annotation in annot.keys():
		if annotation not in annotated:
			notFound.append(annotation)

if len(annot) == 0:
	if args.verbose:
		print("  No annotations were found")
elif len(annotated) == len(annot):
	if args.verbose:
		print("  All '", len(annotated), "' annotations were successfully transferred", sep="")
else:
	print("")
	print("  Warning!! Not all annotations were transferred successfully")
	print("  Annotated nodes in input tree: ", len(annot))
	print("  Annotated nodes in output tree:", len(annotated))
	if len(notFound) > 0 and len(notFound) < 51:
		print("    The following annotations were not found in the tree to be annotated:")
		for a in notFound:
			print("    -", a)
	else:
		tmp = re.sub("\\.[^\\.]+$", "_annotationNotFound.txt", args.treeToAnnot)
		with open(tmp, "w") as tmp1:
			for a in notFound:
				print(a, file=tmp1)
		print("    More than 50 annotations were not found in the tree to be annotated")
		print("    Please check '", tmp,"' for more details", sep="")
	print("")

# Writing files ------------------------------------------------------------------------------------
if args.verbose:
	print("  Exporting annotated tree to:", outFile)
Phylo.write(A, outFile, args.formatOutput)

if args.verbose:
	print("Done")
