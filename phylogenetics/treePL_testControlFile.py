#!/usr/bin/env python3

import argparse
from Bio import Phylo
import re

parser = argparse.ArgumentParser(description="Given a control file, test if every taxon is in the tree and if all calibrations are compatible to one another.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-c", "--controlFile", dest="controlFile", required=True,
					help="A control file for treePL.")

parser.add_argument("-t", "--tree", dest="tree", required=False, default=None,
						  help="A tree file. If not given, will read the tree given in the control file.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format: accepted formats are newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-o", "--out", dest="out", required=False, default=None,
					help="The name of the annotated tree. By default, will add '_calibrated' and the original extension to the input tree file.")

parser.add_argument("-p", "--prune", dest="prune", required=False, action="store_true",
					help="If selected, will prune the tree from taxa not used in the control file.")

parser.add_argument("-k", "--crosscheck", dest="crosscheck", required=False, action="store_true",
					help="If selected, will also check whether the minimum ages are compatible with the maximum ages and viceversa.")

parser.add_argument("-e", "--export", dest="export", required=False, action="store_false",
					help="If selected, will not export an annotated tree.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Reading files ------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading control file")
taxa = {}
taxaList = list()
agesMin = {}
agesMax = {}
for line in open(args.controlFile):
	if args.tree is None:
		if line.startswith("treefile"):
			tree = re.sub(".* ", "", line)
			tree = re.sub("\n", "", tree)
			if args.verbose:
				print("    Tree file detected:", tree)
	if line.startswith("mrca"):
		line = line.strip(" ").split()
		name = str(line[2])
		taxa[name] = list()
		for taxon in range(3, len(line)):
			taxa[name].append(line[taxon])
			taxaList.append(line[taxon])
	elif line.startswith("min"):
		line = line.strip(" ").split()
		agesMin[str(line[2])] = str(line[3])
	elif line.startswith("max"):
		line = line.strip(" ").split()
		agesMax[str(line[2])] = str(line[3])

ages = {}
for key in taxa.keys():
	if key in agesMax and key in agesMin:
		ages[key] = str(agesMax[key] + "-" + agesMin[key])
	elif key in agesMax and key not in agesMin:
		ages[key] = str(agesMax[key] + "-NA")
	elif key not in agesMax and key in agesMin:
		ages[key] = str("NA-" + agesMin[key])
	else:
		print("  ", key, "has no calibrations given")

if args.verbose:
	print("  Reading tree file")
if args.tree is None:
	T = Phylo.read(tree, args.formaTree)
else:
	T = Phylo.read(args.tree, args.formaTree)

if args.out is None:
	if args.tree is None:
		out = re.sub("\\.[^\\.]+$", "_calibrated", tree) + re.sub(".*\\.", ".", tree)
	else:
		out = re.sub("\\.[^\\.]+$", "_calibrated", args.tree) + re.sub(".*\\.", ".", args.tree)
else:
	out = args.export

# Checking nodes -----------------------------------------------------------------------------------
# Check if all taxa in the control file are found in the tree
print("  Checking taxa")
terminals = list()
notInTree = list()
for tip in T.get_terminals():
	terminals.append(tip.name)
for taxon in taxaList:
	if taxon not in terminals:
		notInTree.append(taxon)
if len(notInTree) > 0:
	print("    Warning! The following sequences are not found in the tree:")
	for i in notInTree:
		print("   ", i)
	print("Exiting...")
	import sys
	sys.exit(1)
else:
	print("    All taxa used in calibrations are found in the tree")

# Writing tree -------------------------------------------------------------------------------------
if args.verbose:
	print("  Annotating calibrations in tree")
for key, value in taxa.items():
	clade = T.common_ancestor(value)
	clade.comment = None
	clade.comment = str('[&!color=#FFB000,!name="' + key + "_" + ages[key] +'"]')

if args.export:
	if args.verbose:
		print("  Colouring terminal branches used for calibrations")
	for tip in T.get_terminals():
		tip.comment = None
		if tip.name in taxaList:
			tip.comment = str("[&!color=#FFB000]") # Orange

# check if child nodes are not older than parent nodes
if args.verbose:
	print("  Checking dates")
for clade in T.get_nonterminals():
	if clade.comment is not None:
		nameClade = re.sub('.*name="', "", str(clade.comment))
		nameClade = re.sub('_.*', "", nameClade)
		if nameClade in agesMin:
			ageMin = float(agesMin[nameClade])
		else:
			ageMin = None
		if nameClade in agesMax:
			ageMax = float(agesMax[nameClade])
		else:
			ageMax = None
		for subclade in T.get_path(clade):
			if subclade.comment is not None and subclade.comment != clade.comment:
				nameSubclade = re.sub('.*name="', "", str(subclade.comment))
				nameSubclade = re.sub('_.*', "", nameSubclade)
				if nameSubclade in agesMin and ageMin is not None:
					ageMinSubclade = float(agesMin[nameSubclade])
					if ageMinSubclade < ageMin:
						print("    Warning! ", nameClade, " min age (", ageMin, ") conflicts ", nameSubclade, " min age (",ageMinSubclade, ")", sep="")
						print("      Calibrations: ", nameClade, "(", ages[nameClade],") - ", nameSubclade, "(", ages[nameSubclade],")", sep="")
				if nameSubclade in agesMax and ageMax is not None:
					ageMaxSubclade = float(agesMax[nameSubclade])
					if ageMaxSubclade < ageMax:
						print("    Warning! ", nameClade, " max age (", ageMax, ") conflicts ", nameSubclade, " max age (",ageMaxSubclade, ")", sep="")
						print("      Calibrations: ", nameClade, "(", ages[nameClade],") - ", nameSubclade, "(", ages[nameSubclade],")", sep="")
				if args.crosscheck:
					if nameSubclade in agesMin and ageMax is not None:
						ageMinSubclade = float(agesMin[nameSubclade])
						if ageMinSubclade < ageMax:
							print("    ", nameClade, " max age (", ageMax, ") may conflict ", nameSubclade, " min age (",ageMinSubclade, ")", sep="")
							print("      Calibrations: ", nameClade, "(", ages[nameClade],") - ", nameSubclade, "(", ages[nameSubclade],")", sep="")
					if nameSubclade in agesMax and ageMin is not None:
						ageMaxSubclade = float(agesMax[nameSubclade])
						if ageMaxSubclade < ageMin:
							print("    ", nameClade, " max age (", ageMin, ") may conflict ", nameSubclade, " min age (",ageMaxSubclade, ")", sep="")
							print("      Calibrations: ", nameClade, "(", ages[nameClade],") - ", nameSubclade, "(", ages[nameSubclade],")", sep="")

# prune the tree -----------------------------------------------------------------------------------
if args.prune:
	if args.verbose:
		print("  Prunning", end="")
	toPrune = set()
	for tip in T.get_terminals():
		if tip.comment is None:
			toPrune.add(tip.name)
	toPrunel = len(toPrune)
	i = 0
	pl = 0
	for tip in toPrune:
		if args.verbose:
			i += 1
			p = round(i/toPrunel*100)
			if p > pl:
				pl = p
				print("\r  Pruning ", p, "%",sep="", end="")
		T.prune(tip)

# Export tree ---------------------------------------------------------------------------------------
if args.export:
	import subprocess
	if args.verbose:
		print("  Writing annotated tree to", out)
	Phylo.write(T, out, "nexus")
	subprocess.call(["sed", "-i", "-e",  's/\\\]//g', out])
	subprocess.call(["sed", "-i", "-e",  's/\\[\\\//g', out])

if args.verbose:
	print("Done")
