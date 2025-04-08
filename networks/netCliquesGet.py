#!/usr/bin/env python3.12

import argparse
import networkx as nx
import re

parser = argparse.ArgumentParser(description="Extracts all maximal cliques in a network.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
					help="Net file. A file with three columns, origin node, destination node and identity value.")

parser.add_argument("-n", "--nodes", dest="nodes", required=False, default=None,
					help="A file with the nodes to be searched, each one in one line. This will only yield maximal cliques containing all nodes provided, otherwise it will rise an error.")

parser.add_argument("-s", "--size", dest="size", required=False, default=2, type=int,
					help="A minimal number of nodes for a clique to be considered. By default=2.")

parser.add_argument("-p", "--pattern", dest="pattern", required=False, nargs="+", default=None,
					help="A pattern to be search within the nodes of the clique. At least one node should contain the pattern.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
					help="Output file. By default will add '_subset.net' to input file.")

args = parser.parse_args()

# Setting variables --------------------------------------------------------------------------------
if args.file_out is None:
	out = re.sub("\\.[^\\.]+$", "_cliques.tsv", args.file_in)
else:
	out = args.file_out

# Reading network  ---------------------------------------------------------------------------------
print("  Reading network", flush=True)
G=nx.read_edgelist(args.file_in, delimiter="\t", data=(("id",float),))

# Reading selected nodes, if provided --------------------------------------------------------------
if args.nodes is not None:
	nodes = set()
	for line in open(args.nodes):
		nodes.add(line)
else:
	nodes = None

# Finding cliques  ---------------------------------------------------------------------------------
if nodes is None:
	# cliques = nx.enumerate_all_cliques(G)
	cliques = nx.find_cliques(G)
else:
	cliques = nx.find_cliques(G, nodes)

print("  Writing cliques to", out, flush=True)
with open(out, "w") as outfile:
	count = 0
	exported = 0
	for clique in cliques:
		count += 1
		export = False
		if len(clique) >= args.size:
			export = True
		if args.pattern is not None:
			for p in args.pattern:
				if any(p in n for n in clique):
					export = True
				else:
					export = False
		if export:
			exported += 1
			for n in clique:
				print(n + "\t" + str(exported), file=outfile)

print("    Total cliques:\t", count, sep="", flush=True)
print("    Total exported cliques:\t", exported, sep="", flush=True)
print("Done", flush=True)
