#!/usr/bin/env python3

import argparse
import networkx as nx
import re

parser = argparse.ArgumentParser(description="Extracts connected components that have given node attributes or given string patterns in the nodes.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
					help="Net file. A file with three columns, origin node, destination node and identity value.")

parser.add_argument("-a", "--attributes", dest="attributes", required=False, default=None,
					help="A file with the nodes and the attributes of each node.")

parser.add_argument("-g", "--group", dest="group", required=False, default=None, nargs="+",
					help="The name of the attribute(s) you want to keep.")

parser.add_argument("-p", "--pattern", dest="pattern", required=True, default=None, nargs="+",
					help="A pattern in the nodes to be kept.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
					help="Output file. By default will add '_subset.net' to input file.")

args = parser.parse_args()

if args.file_out is None:
	out = re.sub("\\.[^\\.]+$", "_subset.net", args.file_in)
else:
	out = args.file_out

if args.group is None and args.pattern is None:
	import sys
	print("Error! Please select either a patterns in the nodes or an attribute through the attributes file.")
	sys.exit(0)

if args.group is not None and args.attributes is None:
	import sys
	print("Error! Please select an attributes file along with the attributes group.")
	sys.exit(0)

if args.group is None and args.attributes is not None:
	import sys
	print("Error! Please select an attribute group along with the attributes file.")
	sys.exit(0)

print("  Reading network", flush=True)
G=nx.read_edgelist(args.file_in, delimiter="\t", data=(("id",float),))

nodes=list(G.nodes())

if args.attributes is not None:
	print("  Reading attributes file", end="", flush=True)
	attr = {}
	for line in open(args.attributes):
		line = line[:-1].split("\t")
		key = line[0]
		key = re.sub(" $", "", key)
		value = line[1]
		value = re.sub("^ ", "", value)
		attr[key] = value
	
	print(" and assigning attributes to network", flush=True)
	attributes = {}
	for node in nodes:
		attributes[node] = attr.get(node)
	
	nx.set_node_attributes(G, attributes, "attributes")

print("  Cleaning network", flush=True)
c = 0
k = 0
n = 0
toRemove = set()
count = nx.number_connected_components(G)
CCs = (G.subgraph(CCs) for CCs in nx.connected_components(G))
for CC in CCs:
	c += 1
	print("\r    Working on CC ", c, "/", count, sep="", end="", flush=True)
	groups = nx.get_node_attributes(CC, "attributes")
	groups = groups.values()
	groups = set(groups)
	remove=True
	if args.group is not None:
		for group in args.group:
			if group in groups:
				remove=False
	if args.pattern is not None:
		for node in CC:
			for i in args.pattern:
				if i in node:
					remove = False
					break
	if remove:
		k += 1
		for node in CC:
			n += 1
			toRemove.add(node)

print("\n  Removing", flush=True)
print("    ", k, " CCs", sep="", flush=True)
print("    ", n, " nodes", sep="", flush=True)

G.remove_nodes_from(toRemove)

print("  Writing network to", out, flush=True)
ni = len(nodes)
no = len(list(G.nodes()))
print("    Nodes in:\t", ni, sep="", flush=True)
print("    Nodes out:\t", no, "\t(", round(no/ni*100,1), "%)", sep="", flush=True)
co = nx.number_connected_components(G)
print("    CCs in:\t", count, sep="", flush=True)
print("    CCs out:\t", co, "\t(", round(co/count*100,1), "%)", sep="", flush=True)

nx.write_edgelist(G, out, delimiter="\t", data=["id"])

print("Done", flush=True)
