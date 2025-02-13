#!/usr/bin/env python3

import argparse
import networkx as nx
import statistics as st
import re

parser = argparse.ArgumentParser(description="Analyzes basic properties of a network (number of nodes, number of edges, connectivity, clustering coefficient and number of Connected Components -CCs-), its CCs independently, if any, and its nodes (degree, betweeness, closeness and eccentricity).")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
                          help="Input file. A network file with three columns: 'qseqid sseqid id'.")

parser.add_argument("-o", "--output", dest="file_out", required=False, nargs="+", default=None,
                    help="Output files, in the order (1) network properties, (2) node properties and if present (3) CC properties. By default will add '_stats.log' and '_statsNode.tsv' to the input file after removing the extension (and '_statsCCs.tsv', if more than one CC).")

args = parser.parse_args()

outNetworkCCs = None
if args.file_out is None:
	out = re.sub("\\.[^\\.]+$", "", args.file_in)
	outNetwork = out + "_stats.log"
	outNodes = out + "_statsNode.tsv"
	outNetworkCCs = out + "_statsCCs.tsv"
else:
	if len(args.file_out) < 2:
		print("Error! Please provide at least 2 output names for the network and the nodes properties")
		import sys
		sys.exit(0)
	outNetwork = args.file_out[0]
	outNodes = args.file_out[1]
	if len(args.file_out) >2:
		outNetworkCCs = args.file_out[2]
	if len(args.file_out) > 3:
		print("  Warning! Only the first three output names will be considered")

print("  Reading network", flush=True)
G=nx.read_edgelist(args.file_in, delimiter="\t", data=(("id",float),))

nodes=list(G.nodes())
edges=list(G.edges())

# Splitting the graph into connected components
#CCs = nx.connected_component_subgraphs(G)  # When using a 'networkx' version below 2.1
CCs = (G.subgraph(CCs) for CCs in nx.connected_components(G))
count = 0
for CC in CCs:
    count = count + 1

print("  Calculating properties of the network")
with open(outNetwork, "w") as outfile:
	print(("Input file: " + "\t" + str(args.file_in)), file=outfile, flush=True)
	# Number of nodes
	print(("Number of nodes" + "\t" + str(len(nodes))), file=outfile, flush=True)
	# Number of edges
	print(("Number of edges" + "\t" + str(len(edges))), file=outfile, flush=True)
	# Average number of neighbours 
	tmp = dict(nx.degree(G))
	tmp = tmp.values()
	tmp = st.mean(tmp)
	print(("Connectivity" + "\t" + str(tmp)), file=outfile, flush=True)
	#print(("Connectivity" + "\t" + str(nx.average_degree_connectivity(G))), file=outfile)
	# Clustering coefficient: E/((N(N-1))/2)
	print(("Clustering coefficient" + "\t" + str(nx.density(G))), file=outfile, flush=True)
	# Connected components
	print(("Connected components" + "\t" + str(count)), file=outfile, flush=True)

if count > 1:
	if outNetworkCCs is None:
		print("Error! Please provide a third name for the connected components file")
		import sys
		sys.exit(0)
	print("  Calculating properties of", count, "connected components (CCs)", flush=True)
	#CCs = nx.connected_component_subgraphs(G)  # When using a 'networkx' version below 2.1
	CCs = (G.subgraph(CCs) for CCs in nx.connected_components(G))
	c = 0
	with open(outNetworkCCs, "w") as outfile:
		print("Connected_component \t Nnodes \t Nedges \t Connectivity \t Density", file=outfile, flush=True)
		for CC in CCs:
			c = c +1
			print("\r    Working on CC ", c, "/", count, sep="", end="", flush=True)
			nnodes = len(list(CC.nodes()))
			nedges = len(list(CC.edges()))
			tmp = dict(nx.degree(CC))
			tmp = tmp.values()
			connec = st.mean(tmp)
			densit = nx.density(CC)
			print((str(c) + "\t" + str(nnodes) + "\t" + str(nedges) + "\t" + str(connec) + "\t" + str(densit)), file=outfile, flush=True)

print("\n  Calculating properties of the nodes", flush=True)
#CCs = nx.connected_component_subgraphs(G)  # When using a 'networkx' version below 2.1
CCs = (G.subgraph(CCs) for CCs in nx.connected_components(G))
with open(outNodes, "w") as outfile:
	print("Connected_component \t Node \t Degree \t Betweenness \t Closeness \t Eccentricity", file=outfile, flush=True)
	c = 0
	for CC in CCs:
		c = c + 1
		if count > 1:
			print("\r    Working on CC ", c, "/", count, sep="", end="", flush=True)
		# Degree: Number of edges related to a node
		degree = nx.degree(CC)
		# Betweenness: Frequency of being in a shortest path
		betweeness = nx.betweenness_centrality(CC)
		# Closeness: Average shortest distance between a node and all other nodes
		closeness = nx.closeness_centrality(CC)
		# Eccentricity: Average longest distance between a node and all other nodes
		eccentricity = nx.eccentricity(CC)
		for node in CC.nodes():
			print((str(c) + "\t" + str(node) + "\t" + str(degree[node]) + "\t" + str(betweeness[node]) + "\t" + str(closeness[node]) + "\t" + str(eccentricity[node])), file=outfile, flush=True)

print("\n  Graph properties exported to:\t", outNetwork)
print("  Node properties exported to:\t", outNodes)
if outNetworkCCs is not None:
	print("  CCs properties exported to:\t", outNetworkCCs)
print("Done", flush=True)
