#!/usr/bin/env python3

import argparse
from Bio import Phylo
import re

parser = argparse.ArgumentParser(description="Given an attribute list of the tree tips, will search for conflicting nodes and print a list of potential intruders. Useful for very large trees (>1000 tips). USE WITH CAUTION!!")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

requiredArgs.add_argument("-a", "--attribute", dest="attribute", required=True,
					help="A tab delimited file containing two columns (1) either the exact list of tips or a pattern to be search in the tip name and (2) the attribute of the tip/pattern. If selected 'supergroup', will search for eukaryotic supergroups.")

parser.add_argument("-o", "--output", dest="out", required=False, action="store",
					help="The output file containing the tip names of potential intruders. If not selected, will add '_intruders.list' to the input tree file.")

parser.add_argument("-m", "--min", dest="minimum", required=False, type=int, default=10,
					help="The minimum threshold to be considered an intruder (in percentage: 0-100). If a monophyletic group contains less than 'm' percent of the total tips for a given group, it could be considered an intruder. Default = 10")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-i", "--intruders", dest="intruders", required=False, action="store_true",
					help="If selected, will colour in orange the tips identified as intruders and blue the correct ones to the selected output name replacing or adding '_intruders.tre' extension in nexus format")

parser.add_argument("-p", "--prune", dest="prune", required=False, action="store_true",
					help="If selected, will prune the tree from the selected intruder tips to the selected output name replacing or adding 'Pruned.tre' extension in newick format.")

parser.add_argument("-n", "--none", dest="none", required=False, action="store_false",
					help="If selected, will not ignore tips without an attribute. However if a clade is composed of intruders AND tips without an attribute, the whole clade will be considered an intruder.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# setting variables --------------------------------------------------------------------------------
if args.out is None:
	out = re.sub("\\.[^\\.]+$", "_intruders.list", args.tree)
else:
	out = args.out

if args.intruders:
	treeMarked = re.sub("\\.[^\\.]+$", "_intruders.tre", out)
if args.prune:
	pruned = re.sub("\\.[^\\.]+$", "Pruned.tre", out)

# Reading files ------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading files")
T = Phylo.read(args.tree, args.formaTree)

if args.attribute == 'supergroup' or args.attribute == 'supergroups':
	attribute={"Amoebozoa":     "Amoebozoa",
		  "Obazoa":             "Obazoa",
		  "Nucletmycea":        "Nucletmycea",
		  "Fungi":              "Nucletmycea",
		  "Holozoa":            "Holozoa",
		  "Opisthokonta":       "Holozoa",
		  "Metazoa":            "Holozoa",
		  "Metamonada":         "Metamonada",
		  "CRuMs":              "CRuMs",
		  "Discoba":            "Discoba",
		  "Hemimastigophora":   "Hemimastigophora",
		  "Ancoracysta":        "Ancoracysta",
		  "Malawimona":         "Malawimona",
		  "Ancyromona":         "Ancyromona",
		  "Cryptista":          "Cryptista",
		  "Cryptophyta":        "Cryptista",
		  "Katablepharidophyta":"Cryptista",
		  "Haptista":           "Haptista",
		  "Centroheliozoa":     "Haptista",
		  "Haptophyta":         "Haptista",
		  "Archaeplastida":     "Archaeplastida",
		  "Chloroplastida":     "Archaeplastida",
		  "Chlorophyta":        "Archaeplastida",
		  "Streptophyta":       "Archaeplastida",
		  "Picozoa":            "Archaeplastida",
		  "Telonemia":          "Telonemia",
		  "Rhizaria":           "Rhizaria",
		  "Stramenopiles":      "Stramenopiles",
		  "Colponemida":        "Alveolata",
		  "Alveolata":          "Alveolata"}
else:
	attribute = {}
	for line in open(args.attribute):
		tmp = line.strip().split()
		attribute[tmp[0]]=tmp[1]

# Assigning ----------------------------------------------------------------------------------------
if args.verbose:
	print("  Assigning attributes to tip names")
	if args.none:
		print("    Ignoring tips without attributes")
	else:
		print("    Tips without attributes will be taken into a count")
attrCount = {}
for tip in T.get_terminals():
	comment = None
	if tip.name in attribute.keys():
		comment = attribute[tip.name]
	else:
		for pattern, attr in attribute.items():
			if pattern in tip.name:
				comment = attr
	if args.none and comment is not None:
		if comment not in attrCount.keys():
			attrCount[comment] = 1
		elif comment in attrCount.keys():
			attrCount[comment] += 1
		tip.comment = comment
	elif not args.none:
		if comment is None:
			comment = "None"
		if comment not in attrCount.keys():
			attrCount[comment] = 1
		elif comment in attrCount.keys():
			attrCount[comment] += 1
		tip.comment = comment

# Assigning internal nodes -------------------------------------------------------------------------
if args.verbose:
	print("  Assigning internal nodes")
for clade in T.get_nonterminals():
	unique = set()
	for tip in clade.get_terminals():
		comment = tip.comment
		unique.add(comment)
	if len(unique) == 1:
		clade.comment = next(iter(unique))
	if args.none:
		if len(unique) == 2 and None in unique:
			unique.remove(None)
			clade.comment = next(iter(unique))

# Identifying intruders ----------------------------------------------------------------------------
if args.verbose:
	print("  Identifying intruders")
intruders = set()
prev = None
for clade in T.get_nonterminals():
	if clade.comment is None:
		attrCountClade = {}
		for tip in clade.get_terminals():
			comment = tip.comment
			if args.none and comment is not None:
				if comment not in attrCountClade.keys():
					attrCountClade[comment] = 1
				elif comment in attrCountClade.keys():
					attrCountClade[comment] += 1
			elif not args.none:
				if comment is None:
					comment = "None"
				if comment not in attrCountClade.keys():
					attrCountClade[comment] = 1
				elif comment in attrCountClade.keys():
					attrCountClade[comment] += 1
		for attr, count in attrCountClade.items():
			test = count / attrCount[attr] * 100
			if test < args.minimum:
				for tip in clade.get_terminals():
					if tip.comment == attr and tip.comment != prev:
						intruders.add(tip.name)
		if len(attrCountClade) > 0:
			prev = max(attrCountClade, key=attrCountClade.get)

# Check if there are monophyletic clades with only intruders and tips with no attribute
if args.verbose:
	print("  Identifying intruders without an attribute")
intrudersNone = set()
for clade in T.get_nonterminals():
	tmpIntr = 0
	tmpNone = 0
	tmp = set()
	unique = set()
	for tip in clade.get_terminals():
		if tip.name in intruders:
			tmpIntr += 1
		if tip.comment is None:
			tmpNone += 1
	if (tmpIntr + tmpNone) == clade.count_terminals() and tmpIntr != 0:
		for tip in clade.get_terminals():
			if tip.comment is None:
				intrudersNone.add(tip.name)

if len(intrudersNone) > 0:
	for i in intrudersNone:
		intruders.add(i)

# Exporting list of intruders ----------------------------------------------------------------------
if args.verbose:
	tips = T.count_terminals()
	tipsi = len(intruders)
	print("  From the total", tips, "tips and the given attributes:")
	print("    ", tipsi, " (", round(tipsi / tips * 100, 2), "%) tips were identified as intruders", sep="")
	if len(intrudersNone) > 0:
		print("    of which", len(intrudersNone), "had no attribute but are closely related to intruders")
	print("  Writing list of intruders to", out)
	if tipsi / tips * 100 > args.minimum:
		print("\n    Warning!\n    The proportion of intruders is too high!\n    Please check carefully the attributes file")
		if not args.none:
			print("    Or consider not selecting the '-n/--none' option")
with open(out, "w") as outlist:
	for tip in intruders:
		print(tip, file=outlist)

# Colouring intruders ------------------------------------------------------------------------------
if args.intruders:
	if args.verbose:
		print("  Colouring and exporting input tree with highlighted intruders to", pruned)
	import subprocess
	treeMarked = re.sub("\\.[^\\.]+$", ".tre", out)
	if args.verbose:
		print("    Colouring terminal nodes")
	for clade in T.get_terminals():
		clade.comment = None
		if clade.name in intruders:
			clade.comment = str("[&!color=#FFB000]")
		else:
			clade.comment = str("[&!color=#648FFF]")
	if args.verbose:
		print("    Colouring internal nodes")
	for clade in T.get_nonterminals():
		clade.comment = None
		unique = set()
		for tip in clade.get_terminals():
			unique.add(tip.comment)
		if len(unique) == 1:
			clade.comment = next(iter(unique))
	if args.verbose:
		print("    Exporting")
	Phylo.write(T, treeMarked, "nexus")
	subprocess.call(["sed", "-i", "-e",  's/\\\]//g', treeMarked])
	subprocess.call(["sed", "-i", "-e",  's/\\[\\\//g', treeMarked])

# Pruning intruders --------------------------------------------------------------------------------
if args.prune:
	if args.verbose:
		print("  Pruning")
	for tip in intruders:
		T.prune(tip)
	if args.verbose:
		print("    Writing pruned tree to", pruned)
	for clade in T.get_terminals():
		clade.comment = None
	for clade in T.get_nonterminals():
		clade.comment = None
	Phylo.write(T, pruned, "newick")

if args.verbose:
	print("Done")
