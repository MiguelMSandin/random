#!/usr/bin/env python3

import argparse
import networkx as nx
import statistics as st
import re

parser = argparse.ArgumentParser(description="Analyzes basic properties of a network (number of nodes, number of edges, connectivity, clustering coefficient and number of connected components), its connected components independently, if any, and its nodes (degree, betweeness, closeness and eccentricity).")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input file. A network file with three columns: 'qseqid sseqid id'.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
                    help="Output file. Returns two files: 'output_network.log' and 'output_nodes.log' with the properties of the network and the statistics for each node respectively. By default will add '_stats.log' and '_nodeStats.log' to the input file after removing the extension.")

args = parser.parse_args()

if args.file_out is None:
	out = re.sub("\\.[^\\.]+$", "", args.file_in)
else:
	out = args.file_out

outNetwork = out + "stats.log"
outNodes = out + "_nodeStats.log"

print("\n  Reading network", flush=True)
G=nx.read_edgelist(args.file_in, delimiter="\t", data=(("id",float),))

nodes=list(G.nodes())
edges=list(G.edges())

# Splitting the graph into connected components
#CCs = nx.connected_component_subgraphs(G)  # When using a 'networkx' version below 2.1
CCs = (G.subgraph(CCs) for CCs in nx.connected_components(G))
count = 0
for CC in CCs:
    count = count + 1

print("\n  Calculating properties of the network")
with open(outNetwork, "w") as outfile:
    print(("Input file: " + "\t" + str(args.file_in)), file=outfile, flush=True)
    # Number of nodes
    print("    Number of nodes", end="", flush=True)
    print(("Number of nodes" + "\t" + str(len(nodes))), file=outfile, flush=True)
    # Number of edges
    print(": done \n    Number of edges", end="", flush=True)
    print(("Number of edges" + "\t" + str(len(edges))), file=outfile, flush=True)
    # Average number of neighbours 
    print(": done \n    Connectivity", end="", flush=True)
    tmp = dict(nx.degree(G))
    tmp = tmp.values()
    tmp = st.mean(tmp)  
    print(("Connectivity" + "\t" + str(tmp)), file=outfile, flush=True)
    #print(("Connectivity" + "\t" + str(nx.average_degree_connectivity(G))), file=outfile)
    # Clustering coefficient: E/((N(N-1))/2)
    print(": done \n    Clustering coefficient", end="", flush=True)
    print(("Clustering coefficient" + "\t" + str(nx.density(G))), file=outfile, flush=True)
    # Connected components
    print(": done \n    Connected components", end="", flush=True)
    print(("Connected components" + "\t" + str(count)), file=outfile, flush=True)
    print(": done\n", flush=True)

if count > 1:
    print("  There are '", count, "' connected components (CC). Analyzing each one of them individually as an independent network:", sep="", flush=True)
    outNetworkCCs = out + "_network_CCs"
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
        print("\n", flush=True)

print("  Calculating properties of the nodes", flush=True)
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
    print("\nDone", flush=True)
