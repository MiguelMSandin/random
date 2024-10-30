#!/usr/bin/env python3

import argparse
import networkx as nx

parser = argparse.ArgumentParser(description="Calculates the shortest path (number of nodes) that requires to go from state 'from' to state 'to' of the given attribute for every node of state 'from'.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
						  help="Input file. A network file with three columns: 'qseqid sseqid id'.")

requiredArgs.add_argument("-a", "--attributeFile", dest="file_attr", required=True,
						  help="Attributes file. A file with the nodes in the first column and the binary attribute in the second column. Columns separated by '\\t'.")

parser.add_argument("-b", "--binaryStates", dest="binary_states", required=False, nargs=2, default=None,
					help="The binary states of the attribute to look for the shortest path from 'from' to 'to' entered as a string. If not provided will take the first appearing state as 'from' and the second state as 'to'.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
					help="Output file. A file with the shortest path of every node labelled as 'from' to a node labelled as 'to'. By default, will add '_shortest.tsv' to the input file (after deleting the extension)")

parser.add_argument("-r", "--header", dest="header", required=False, action="store_false",
					help="If selected, will remove the headers (first row) of the attribute file.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

if args.file_out is None:
	import re
	out = re.sub("\\.[^\\.]+$", "_shortest.tsv", args.file_in)
else:
	out = args.file_out

if args.verbose:
	print("  Reading network file")
G=nx.read_edgelist(args.file_in, delimiter="\t", data=(("id",float),))

if args.verbose:
	print("  Reading attributes file")
with open(args.file_attr) as attributes:
	attr = {}
	i = 0
	if args.binary_states is None:
		states = set()
	for line in attributes:
		i += 1
		if (args.header is False) or (args.header and i != 1):
			line = line.strip().split('\t')
			key = line[0]
			val = line[1]
			attr[key] = val
			if args.binary_states is None:
				states.add(val)

# Setting the names of the states to look from 'from' to 'to'
if args.binary_states is None:
	if len(states) < 2:
		print(" Warning!! There are less than two states in the attributes file. This will yield an error...")
	if len(states) > 2:
		print(" Warning!! There are more than two states in the attributes file. Some states will be ignored...")
	states_from = list(states)[0]
	states_to = list(states)[1]
else:
	states = args.binary_states
	states_from = states[0]
	states_to = states[1]

if args.verbose:
	print("  Setting attributes and states")
# Now matching order of every node with the attribute
group = {}
for node in G.nodes():
	group[node] = attr.get(node)
nx.set_node_attributes(G, group, "groups")
lfrom = len([n for n in G.nodes() if G.nodes[n]["groups"] == states_from])
lto = len([n for n in G.nodes() if G.nodes[n]["groups"] == states_to])

if args.verbose:
	print("  Nodes in network: \t", len(G.nodes()),
	   "\n    Nodes with state 'from' (attribute name '", states_from, "'): \t", lfrom, 
	   "\n    Nodes with state 'to'   (attribute name '", states_to, "'): \t", lto, sep="")

if args.verbose:
	print("  Calculating")
#CCs = nx.connected_component_subgraphs(G)  # When using a 'networkx' version below 2.1
CCs = (G.subgraph(CCs) for CCs in nx.connected_components(G))
with open(out, "w") as outfile:
	print("cc\tnode_from\tshortest\tnode_to", file=outfile)
	c = 0
	i = 0
	for CC in CCs:
		c = c + 1
		nfrom = [n for n in CC.nodes() if CC.nodes[n]["groups"] == states_from]
		nto = [n for n in CC.nodes() if CC.nodes[n]["groups"] == states_to]
		for n in nfrom:
			i = i + 1
			if args.verbose:
				print("\r    Working on node ", i, "/", lfrom, sep="", end="")
			if states_to in nx.get_node_attributes(CC, "groups").values():
				#shortest = sorted([nx.shortest_path(CC, n, ntoi) for ntoi in nto], key=len)[0]
				shortest = sorted([nx.shortest_path(CC, n, ntoi) for ntoi in nto], key=lambda x: len(x))[0]
				end = shortest[len(shortest)-1]
				shortest = len(shortest)-1
			else:
				shortest = "Inf"
				end = "Node_to_not_in_CC"
			print((str(c) + "\t" + str(n) + "\t" + str(shortest) + "\t" + str(end)), file=outfile)

if args.verbose:
	print("\nDone")
