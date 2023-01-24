#!/usr/bin/env python3

import argparse
import networkx as nx
import re

parser = argparse.ArgumentParser(description="Extracts connected components that have given node attributes to speed processing.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Net file. A file with three columns, origin node, destination node and identity value.")

requiredArgs.add_argument("-a", "--attributes", dest="attributes", required=True,
                    help="A file with the nodes and the attributes of each node.")

requiredArgs.add_argument("-g", "--group", dest="group", required=True, nargs="+",
                    help="The name of the attribute(s) you want to keep.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
                    help="Output file. By default will add '_subset.net' to input file.")

args = parser.parse_args()

if args.file_out is None:
	out = re.sub("\\.[^\\.]+$", "_subset.net", args.file_in)
else:
	out = args.out

print("  Reading network")
G=nx.read_edgelist(args.file_in, delimiter="\t", data=(("id",float),))

nodes=list(G.nodes())

print("  Reading attributes file", end="")
attr = {}
for line in open(args.attributes):
    line = line[:-1].split("\t")
    key = line[0]
    key = re.sub(" $", "", key)
    value = line[1]
    value = re.sub("^ ", "", value)
    attr[key] = value

print(" and assigning attributes to network")
attributes = {}
for node in nodes:
    attributes[node] = attr.get(node)

nx.set_node_attributes(G, attributes, "attributes")

print("  Cleaning network")
c = 0
k = 0
n = 0
toRemove = set()
count = nx.number_connected_components(G)
CCs = (G.subgraph(CCs) for CCs in nx.connected_components(G))
for CC in CCs:
    c += 1
    print("\r    Working on CC ", c, "/", count, sep="", end="")
    groups = nx.get_node_attributes(CC, "attributes")
    groups = groups.values()
    groups = set(groups)
    remove=True
    for group in args.group:
        if group in groups:
            remove=False
    if remove:
        k += 1
        for node in CC:
            n += 1
            toRemove.add(node)

print("\n  Removing:")
print("      ", k, " CCs", sep="")
print("      ", n, " nodes", sep="")

G.remove_nodes_from(toRemove)

print("    Initial network had  ", len(nodes), " nodes and ", count, " CCs", sep="")
print("    Exported network has ", len(list(G.nodes())), " nodes and ", nx.number_connected_components(G), " CCs", sep="")

print("  Writing network to", out)
nx.write_edgelist(G, out, delimiter="\t", data=["id"])

print("Done")
