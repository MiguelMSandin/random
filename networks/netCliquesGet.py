#!/usr/bin/env python3.12

import argparse
import networkx as nx
import re

parser = argparse.ArgumentParser(description="Extracts all maximal cliques in a network, and exports a tab delimited table with the node and the given clique they belong to, ordered arbitrarily.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
					help="Net file. A file with three columns, origin node, destination node and identity value.")

parser.add_argument("-s", "--size", dest="size", required=False, default=3, type=int,
					help="A minimal number of nodes for a clique to be considered. By default=3.")

parser.add_argument("-p", "--pattern", dest="pattern", required=False, nargs="+", default=None,
					help="A pattern to be search within the nodes of the clique. At least one node should contain the pattern.")

parser.add_argument("-a", "--all", dest="alln", required=False, default=None, action="store_true",
					help="This option will force that all nodes contain the given pattern(s), if provided.")

parser.add_argument("-n", "--nodes", dest="nodes", required=False, default=None,
					help="A file with the nodes to be searched, each one in one line. This will only search cliques within the given nodes, ignoring the '-p/--pattern' argument if provided. This will considerably decrease computing time in large networks.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
					help="Output file. By default will add '_subset.net' to input file.")

args = parser.parse_args()

# Setting variables --------------------------------------------------------------------------------
if args.file_out is None:
	out = re.sub("\\.[^\\.]+$", "_cliques.tsv", args.file_in)
else:
	out = args.file_out

print("    Minimal clique size:", args.size, flush=True)
if args.pattern is not None:
	print("    Cliques containing: ", *args.pattern, flush=True)
	if args.alln is not None:
		print("      All nodes should contain these patterns", flush=True)

# Reading network  ---------------------------------------------------------------------------------
print("  Reading network", flush=True)
G=nx.read_edgelist(args.file_in, delimiter="\t", data=(("id",float),))

# Reading selected nodes, if provided --------------------------------------------------------------
if args.nodes is not None:
	print("  Reading selected nodes", flush=True)
	import itertools
	nodes = set()
	removed = 0
	for line in open(args.nodes):
		line = line.rstrip()
		if line in G.nodes():
			nodes.add(line)
		else:
			removed += 1
	nodes = list(nodes)
	if removed > 0:
		print("   ", removed, "nodes were not found in the graph, and therefore ignored")
	lnodes = len(nodes)
	print("    Nodes considered:   ", lnodes)
else:
	nodes = None

# Finding cliques  ---------------------------------------------------------------------------------
print("  Writing cliques to:", out, flush=True)
with open(out, "w") as outfile:
	exported = 0
	if nodes is None:
		count = 0
		cliques = nx.find_cliques(G)
		for clique in cliques:
			count += 1
			if len(clique) >= args.size:
				testSize = True
			else:
				testSize = False
			if args.pattern is not None:
				for p in args.pattern:
					if args.alln:
						if all(p in n for n in clique):
							testPattern = True
						else:
							testPattern = False
					else:
						if any(p in n for n in clique):
							testPattern = True
						else:
							testPattern = False
			else:
				testPattern = True
			if testSize and testPattern:
				exported += 1
				for n in clique:
					print(n + "\t" + str(exported), file=outfile)
		print("    Total cliques found:\t", count, sep="", flush=True)
	else:
		print("    Finding neighbours", flush=True)
		neighbours = list()
		for node in nodes:
			new = True
			tmp = list()
			tmp.append(node)
			for t in list(G.neighbors(node)):
				tmp.append(t)
			tmp.sort()
			for n in neighbours:
				if n == tmp:
					new = False
			if new:
				neighbours.append(tmp)
		neighbours.sort(key=lambda x:[len(x), x[0]], reverse=True)
		exported = 0
		ni = 0
		for n in neighbours:
			ni += 1
			tmp = round(ni/len(neighbours)*100)
			print("\r    Checking for cliques\t", tmp, "%", sep="", end="", flush=True)
			try:
				clique = list(nx.find_cliques(G, list(n)))
			except ValueError:
				clique = None
			if clique is not None:
				clique = list(itertools.chain(*clique))
				export = True
				for node in clique:
					if node not in nodes:
						export = False
				if export:
					exported += 1
					for node in clique:
						print(node + "\t" + str(exported), file=outfile, flush=True)
		print("", flush=True)

print("    Total exported cliques:\t", exported, sep="", flush=True)
print("Done", flush=True)
