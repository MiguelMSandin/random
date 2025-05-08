#!/usr/bin/env python3

import argparse
import networkx as nx
import re

parser = argparse.ArgumentParser(description="Extracts louvain clusters from within a graph. For more details see 'louvain_communities' from the 'networkx' module.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
                          help="A network input file.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
                    help="Output filename. By default will replace the extension by '_louvain.net'.")

parser.add_argument("-r", "--resolution", dest="resolution", required=False, default=None,
					help="The resolution of the clusters. By default = 1. If r < 1, the algorithm favors larger communities, if r > 1 favors smaller communities.")

args = parser.parse_args()

outNetworkCCs = None
if args.file_out is None:
	outFile = re.sub("\\.[^\\.]+$", "_louvain.net", args.file_in)
else:
	outFile = args.file_out

print("  Reading network", flush=True)
G=nx.read_edgelist(args.file_in, delimiter="\t", data=False)

nodes=list(G.nodes())
edges=list(G.edges())

print("  Calculating louvain communities", flush=True)
louvain=nx.community.louvain_communities(G)

print("  Exporting louvain communities to:", outFile, flush=True)
with open(outFile, "w") as outfile:
	for communityi in louvain:
		tmp = set()
		for node in communityi:
			neighbours = G.neighbors(node)
			for n in neighbours:
				if n in communityi and n != node:
					out = node + "\t" + n
					outr = n + "\t" + node
					if out not in tmp and outr not in tmp:
						tmp.add(out)
						print(out, file=outfile)

print("Done", flush=True)
