#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Phylo
import re
import subprocess

parser = argparse.ArgumentParser(description="Colours a tree based on a table with colours and exports a coloured nexus tree file.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

requiredArgs.add_argument("-c", "--colours", dest="colours", required=True,
					help="A tab delimited file containing two columns (1) either the exact list of tips to be coloured or a pattern to be search in the tip name and (2) the colour of the tip/pattern. If selected 'eukProt', will automatically generate eukProt colours.")

parser.add_argument("-o", "--output", dest="out", required=False, action="store",
					help="The output tree file. If not selected, will add '_coloured' to the input tree file.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

#parser.add_argument("-a", "--action", dest="action", required=False, default='branches',
					#choices=['b', 'branches', 't', 'tips', 'a', 'all'],
					#help="Either colouring only the branches (default: branches), the tip names (tips), or both (all). The initials also work (b, t and a respectively). NOT YET IMPLEMENTED.")

parser.add_argument("-i", "--internal", dest="internal", required=False, action="store_false",
					help="If selected, will not colour internal branches.")

parser.add_argument("-n", "--none", dest="none", required=False, action="store_false",
					help="If selected, will not ignore tips without an attribute. This is useful when colouring internal nodes, so internal branches will not be coloured if one child has no attribute.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# setting variables --------------------------------------------------------------------------------
if args.out is None:
	out = re.sub("\\.[^\\.]+$", "_coloured.", args.tree) + re.sub(".*\\.", "", args.tree)
else:
	out = args.out

# Reading files ------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading files")

T = Phylo.read(args.tree, args.formaTree)

if args.colours == 'eukProt' or args.colours == 'EukProt' or args.colours == 'eukprot':
	colours={"Amoebozoa":       "#9ecae1",
		  "Breviatea":          "#6baed6",
		  "Apusomonadida":      "#6baed6",
		  "Obazoa":             "#6baed6",
		  "Nucletmycea":        "#4292c6",
		  "Rotosphaerida":      "#4292c6",
		  "Fungi":              "#4292c6",
		  "Holozoa":            "#2171b5",
		  "Ichthyosporea":      "#2171b5",
		  "Choanoflagellata":   "#2171b5",
		  "Metazoa":            "#2171b5",
		  "Metamonada":         "#41ab5d",
		  "CRuMs":              "#8dd3c7",
		  "Discoba":            "#fdae6b",
		  "Hemimastigophora":   "#dadaeb",
		  "Ancoracysta":        "#dadaeb",
		  "Malawimona":         "#bdbdbd",
		  "Ancyromona":         "#969696",
		  "Cryptista":          "#fa9fb5",
		  "Cryptophyta":        "#fa9fb5",
		  "Katablepharidophyta":"#fa9fb5",
		  "Haptista":           "#ffffb3",
		  "Centroheliozoa":     "#ffffb3",
		  "Haptophyta":         "#ffffb3",
		  "Archaeplastida":     "#ccebc5",
		  "Picozoa":            "#ccebc5",
		  "Telonemia":          "#bcbddc",
		  "Rhizaria":           "#9e9ac8",
		  "Stramenopiles":      "#807dba",
		  "Colponemida":        "#6a51a3",
		  "Alveolata":          "#6a51a3"}
else:
	colours = {}
	for line in open(args.colours):
		tmp = line.strip().split()
		colours[tmp[0]]=tmp[1]

# Colouring ----------------------------------------------------------------------------------------
if args.verbose and args.internal:
	print("  Colouring terminal nodes")
if args.verbose and not args.internal:
	print("  Colouring")

c = 0
C = 0
for tip in T.get_terminals():
	coloured = False
	if tip.name in colours.keys():
		colour = colours[tip.name]
		if "#" not in colour:
			colour = "#" + colour
		tip.comment = str("[&!color=" + colour + "]")
		coloured = True
	else:
		for pattern, colour in colours.items():
			if pattern in tip.name:
				if "#" not in colour:
					colour = "#" + colour
				tip.comment = str("[&!color=" + colour + "]")
				coloured = True
	if coloured:
		c += 1
	else:
		C += 1

# Colouring internal nodes -------------------------------------------------------------------------
i = 0
I = 0
if args.internal:
	if args.verbose:
		print("  Colouring internal nodes, this might take a while...")
	for clade in T.get_nonterminals():
		coloured = False
		unique = set()
		for tip in clade.get_terminals():
			comment = tip.comment
			unique.add(comment)
		if len(unique) == 1:
			clade.comment = next(iter(unique))
			coloured = True
		if args.none:
			if len(unique) == 2 and None in unique:
				unique.remove(None)
				clade.comment = next(iter(unique))
				coloured = True
		if coloured:
			i += 1
		else:
			I += 1

# Writing file -------------------------------------------------------------------------------------
if args.verbose:
	if c > 0:
		print("  In total", str(i+c), "branches were coloured:")
		print("    Of which ", c, " are terminal branches (", round(c/(c+C)*100,2), "% coloured)", sep="")
		if args.internal:
			print("    Of which ", i, " are internal branches (", round(i/(i+I)*100,2), "% coloured)", sep="")
	else:
		print("    0 branches were coloured, please check table for possible typos.")
	print("  Writting file to:", out)
Phylo.write(T, out, "nexus")
subprocess.call(["sed", "-i", "-e",  's/\\\]//g', out])
subprocess.call(["sed", "-i", "-e",  's/\\[\\\//g', out])

if args.verbose:
	print("Done")
