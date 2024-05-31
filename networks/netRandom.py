#!/usr/bin/env python3

import argparse
import random

parser = argparse.ArgumentParser(description="Generates a random network. See the arguments for further details:",
formatter_class=argparse.RawDescriptionHelpFormatter,
epilog='''\
Examples:
A random network with 20 nodes and 100 edges exported to random.net:
$ netRandom.py

A random network with 200 nodes and 10000 edges exported to random_network.tsv:
$ netRandom.py -n 200 -e 10000 -o random_network.tsv

A random network with 20 nodes and 100 edges, a 0.5 targetted assortativity, 3 states named 'A', 'B' and 'C' and with 0.5, 0.25 and 0.25 occurrence proportions respectively for each state:
$ netRandom.py -a 0.5 -s A B C -p 0.5 0.25 0.25
''')

# Add the arguments to the parser

parser.add_argument("-n", "--nodes", dest="node_number", required=False, default=20,
					help="The number of nodes. By default=20")

parser.add_argument("-e", "--edges", dest="edge_number", required=False, default=100,
					help="The number of edges. By default=100")

parser.add_argument("-a", "--assortativity", dest="assortativity", required=False, type=float, default=0,
					help="The targeted assortativity: A measure of the preference for nodes with a given attribute in a network to attach to other nodes with identical attributes. Default=0, only values [-1 ,0] are considered.")

parser.add_argument("-s", "--states", dest="states", required=False, nargs="+", default=["A", "B"],
					help='The state names. Only used if assortativity different than 0. The attribute will be used for the node names followed by increasing integers. By default=["A", "B"].')

parser.add_argument("-p", "--proportions", dest="proportions", required=False, nargs="+", type=float, default=[1],
					help='The proportion of the states. Only used if assortativity different than 0. If only one proportion is given, it will assume that proportion to all states. By default=1')

parser.add_argument("-c", "--headers", dest="headers", required=False, default=None,
					help="The column names. By default will not print column names.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default="random.net",
					help="The output file.")

parser.add_argument("-i", "--information", dest="information", required=False, action="store_true",
                    help="If selected, will export a table with the nodes and its attribute to 'output`_info.tsv, after removing the extension, if any, and will also print to the console basic information of the randomization process and the final assortativity. Only used if assortativity different than 0.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
                    help="If selected, will not print information in the console.")

args = parser.parse_args()

# Start generating the random net
if args.assortativity == 0:
	if args.verbose:
		print("  Exporting a net of", args.node_number, "nodes and", args.edge_number, "edges to", args.file_out)
	nodes = []
	for i in range(args.node_number):
		tmp = str("n" + str(i+1))
		nodes.append(tmp)
	with open(args.file_out, "w") as outfile:
		for i in range(args.edge_number):
			tmp = random.sample(nodes, 2)
			print(tmp[0] + "\t" + tmp[1], file=outfile)
else:
	# Check parsed parameters
	if args.assortativity > 1:
		assortativity = 1
		print("  Warning! Assortativity can only be [-1, 1]. I have changed the value to '1'")
	elif args.assortativity < -1:
		assortativity = -1
		print("  Warning! Assortativity can only be [-1, 1]. I have changed the value to '-1'")
	else:
		assortativity = args.assortativity
	if len(args.states) != len(args.proportions) & len(args.proportions) != 1:
		print("  Error: The number of states (", len(args.states), ") is different than the number of proportions (", len(args.proportions), ") provided\nExiting", sep="")
		import sys
		sys.exit(1)
	if args.verbose:
		print("  Exporting a net of", args.node_number, "nodes,", args.edge_number, "edges and", assortativity, "targetted assortativity to", args.file_out)
	node_number = int(args.node_number)
	edge_number = int(args.edge_number)
	#
	# First get the states
	states = args.states
	# Now get the proportions
	proportions = []
	if len(args.proportions) == 1:
		for i in states:
			proportions.append(args.proportions[0])
	else:
		proportions = args.proportions
	#
	# Generate a random dictionary of nodes for the given attributes and proportions
	statesp = {}
	for i in states:
		statesp[i] = []
	for i in range(node_number):
		attribute = random.choices(states, weights=proportions, k=1)
		attribute = attribute[0]
		n = len(statesp[attribute])
		statesp[attribute].append(str(attribute + str(n+1)))
	#
	# Now start creating the random network
	if args.information:
		import networkx as nx
		net = nx.Graph()
	with open(args.file_out, "w") as outfile:
		for i in range(edge_number):
			# Generate a random attribute
			attribute = random.choices(states, weights=proportions, k=1)
			attribute = attribute[0]
			probs = []
			for i in states:
				if i == attribute:
					tmp = (assortativity+1)/2
				else:
					tmp = 1-((assortativity+1)/2)
				probs.append(tmp)
			attribute2 = random.choices(states, weights=probs, k=1)
			attribute2 = attribute2[0]
			if attribute == attribute2:
				nodes = random.choices(statesp[attribute], k=2)
			else:
				nodes = random.choices(statesp[attribute], k=1)
				tmp = random.choices(statesp[attribute2], k=1)
				nodes.append(tmp[0])
			print(nodes[0] + "\t" + nodes[1], file=outfile)
			if args.information:
				net.add_edge(nodes[0], nodes[1])

if args.information & args.assortativity != 0:
	import re
	tmp = re.sub("\\.[^\\.]+$", "_info.tsv", args.file_out)
	if args.verbose:
		print("  Exporting attribute table to", tmp)
	attr = {}
	with open(tmp, "w") as outfile:
		for k in statesp.keys():
			tmp = statesp[k]
			for v in tmp:
				print(v + "\t" + k, file=outfile)
				attr[v] = k
	print("  After randomization there are:")
	for k,v in statesp.items():
		tmp = len(v)
		print("    -", tmp, " nodes of the attribute ", k, " (", round(tmp/node_number, 2), "%)", sep="")
	nx.set_node_attributes(net, attr, "groups")
	final_assortativity = nx.attribute_assortativity_coefficient(net, "groups")
	print("    Final assortativity:", final_assortativity)

if args.verbose:
	print("Done")
