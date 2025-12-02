#!/usr/bin/env python3

import argparse
from Bio import Phylo

parser = argparse.ArgumentParser(description="Given a tree and a list of names be annotated, it will find those names in the tip names and annotate it in the last common ancestor of the clade.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-t", "--tree", dest="tree", required=True,
						  help="A tree file.")

parser.add_argument("-l", "--list", dest="listNames", required=False, default=None,
					help="A file with each name in one line.")

parser.add_argument("-n", "--names", dest="names", required=False, default=None, nargs="+",
					help="All names separated with a space.")

parser.add_argument("-f", "--format", dest="formaTree", required=False, default='newick',
					help="The tree file format. Accepted formats are: newick (default), nexus, nexml, phyloxml or cdao.")

parser.add_argument("-o", "--output", dest="output", required=False, default=None,
					help="The output name of the annotated tree file. By default it will add '_annot' before the extension of the input tree file.")

parser.add_argument("-F", "--formatOut", dest="formatout", required=False, default='nexus',
					help="The annotated tree file format. Accepted formats are: newick, nexus (default), nexml, phyloxml or cdao.")

parser.add_argument("-e", "--error", dest="error", required=False, type=float, default=0,
					help="Allows an error in the tip names, so that not all tip names have to contained the name and the annotated clade do not have to be the last common ancestor of all named tips. Default=%(default)s, it will take the last common ancestor node. If, e.g., 0.95 is chosen means that the first node containing at least 95%% of the tips with the given annotations it will be annotated.")

parser.add_argument("-c", "--collapse", dest="collapse", required=False, action="store_true",
					help="If selected, will collapse annotated clades. Warning! This option is only tested with nexus as an output format and on non-ultrametric trees.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Setting variables --------------------------------------------------------------------------------
# Output file
if args.output is None:
	import re
	outFile = re.sub("\\.[^\\.]+$", "_annot.", args.tree) + re.sub(".*\\.", "", args.tree)
else:
	outFile = args.output

if args.collapse:
	import statistics

# Reading files -------------------------------------------------------------------------------------
if args.listNames is None and args.names is None:
	print("  Please provide annotations through '-l/--list' and/or '-n/--names' option.")
	import sys
	sys.exit(0)
else:
	annotations = {}
	w = 0
	if args.listNames is not None:
		for line in open(args.listNames):
			annotations[line]=False
	if args.names is not None:
		for line in args.names:
			if line in annotations.keys():
				w += 1
			else:
				annotations[line]=False
	if w > 0:
		print("  Warning,", str(w),"nnotations were provided by a list and by names and therefore they were dereplicated.")

if args.verbose:
	print("  Reading tree:  ", args.tree)
tree = Phylo.read(args.tree, args.formaTree)

# Loop through the tree and annotate ---------------------------------------------------------------
if args.error > 0:
	i = 0
	for annot in annotations.keys():
		i += 1
		print("\r    Annotating ", annot, " (", i, "/", str(len(annotations)), ")", sep="", end="\n")
		for clade in tree.get_nonterminals():
			if not annotations[annot]:
				counta=0
				for tip in clade.get_terminals():
					if annot in tip.name:
						counta += 1
				prop = counta / clade.count_terminals()
				if prop >= args.error:
					if args.collapse:
						bl = list()
						for t in clade.get_terminals():
							bl.append(clade.distance(t))
						bl = str(round(statistics.mean(bl), 2))
						clade.comment = str('&!name="' + annot + '",!collapse={"collapsed",' + bl + "}")
						annotations[annot] = True
					else:
						clade.comment = str('&!name="' + annot + '"')
						annotations[annot] = True
					continue
			continue
	for annot, val in annotations.items():
		if not val:
			print("      Annotation", annot, "not found")
else:
	print("  Extracting Last Common Ancestor nodes")
	lca = {}
	for tip in tree.get_terminals():
		for annot in annotations.keys():
			if annot in tip.name:
				if annot not in lca.keys():
					lca[annot] = list()
					lca[annot].append(tip.name)
				else:
					lca[annot].append(tip.name)
	print("  Annotating")
	i=0
	for annot,val in lca.items():
		i += 1
		print("\r    Annotating ", annot, " (", i, "/", str(len(annotations)), ")", sep="", end="\n")
		clade = tree.common_ancestor(val)
		if args.collapse:
			bl = list()
			for t in clade.get_terminals():
				bl.append(clade.distance(t))
			bl = str(round(statistics.mean(bl), 2))
			clade.comment = str('&!name="' + annot + '",!collapse={"collapsed",' + bl + "}")
		else:
			clade.comment = str('&!name="' + annot + '"')

# Writing files ------------------------------------------------------------------------------------
if args.verbose:
	print("  Exporting tree to:", outFile)
Phylo.write(tree, outFile, args.formatout)

if args.verbose:
	print("Done")
