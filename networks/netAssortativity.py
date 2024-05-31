#!/usr/bin/env python3

import argparse
import networkx as nx
import statistics as st
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Analyzes the assortativity of a network by groups.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
					help="Input file. A network file with three columns: 'qseqid sseqid id'.")

requiredArgs.add_argument("-a", "--attributes", dest="file_attr", required=True,
					help="Attributes file. A file with the nodes in the first column and the rest of the groups in the following columns matching the given node, separated by '\\t'. Note: Remember to add the column names, so it takes that name to sort out the output.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
					help="Output file. Returns a file with the output. If there are more than 1 connected components (CC), another file will be exported to 'file_out_CCs' with the assortativity calculated independently for every CC. By default, will add '_assortativity'")

parser.add_argument("-c", "--headers", dest="headers", required=False, action="store_false",
					help="If selected, will skip the first row of the network interpreted as headers.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

if args.file_out is None:
        out = re.sub("\\.[^\\.]+$", "_assortativity", args.file_in)
else:
        out = args.file_out

if args.verbose:
	print("  Reading network")
if args.headers:
	G=nx.read_edgelist(args.file_in, delimiter="\t", data=(("id",float),))
else:
	with open(args.file_in, 'rb') as net:
		next(net, '')
		G = nx.read_edgelist(net, delimiter='\t', data=(("id",float),))

nodes=list(G.nodes())
edges=list(G.edges())

if args.verbose:
	print("  Reading attributes file")

attributes = pd.read_csv(args.file_attr, sep="\t")

# Opening files
with open(out, "w") as outfile:
	print("group \t assortativity_groups \t attributes", file=outfile)
# Check if there are more than one connected commponents (CCs)
CCs = (G.subgraph(CCs) for CCs in nx.connected_components(G)) # When using a 'networkx' version above 2.4
#CCs = nx.connected_component_subgraphs(G)   # When using a 'networkx' version below 2.1
count = 0
for CC in CCs:
	count = count + 1
if count > 1:
	outNetworkCCs = out + "_CCs"
	with open(outNetworkCCs, "w") as outfile:	
		print("group \t connected_component \t assortativity_groups \t attributes", file=outfile)

for i in list(range(1, len(attributes.columns))):
	tmp = pd.concat([attributes.iloc[:,0], attributes.iloc[:,i]], axis=1)
	if args.verbose:
		print("  Working on attribute '", list(tmp.columns)[1], "' (", i, "/", len(attributes.columns)-1, ")", sep="")
	# Creating a dictionary from the attributes file 
	attr = {}
	for line in tmp.index:
		key = tmp.iloc[:,0][line]
		val = tmp.iloc[:,1][line]
		if str(val) == "nan":
			val = "NA"
		attr[key] = val
	
	# Now matching order of every node with the attribute
	group = {}
	for node in nodes:
		group[node] = attr.get(node)
	
	# And adding the attributes to the network
	nx.set_node_attributes(G, group, "groups")
	
	if args.verbose:
		print("    Calculating: assortativity by groups", end="")
	assortAttr = nx.attribute_assortativity_coefficient(G, "groups")
	if str(assortAttr) == "nan":
		assortAttr = 0
	if args.verbose:
		print(" and unique groups:", end="")
	groups = nx.get_node_attributes(G, "groups")
	groups = groups.values()
	groups = set(groups)
	if None in groups:
		groups = set(list(filter(None, groups)))
	groups = sorted(groups)
	groups = "|".join(groups)
	
	with open(out, "a") as outfile:
		print((str(list(tmp.columns)[1]) + "\t" + str(assortAttr) + "\t" + str(groups)), file=outfile)
	if args.verbose:
		print(" exported")
	
	if count > 1:
		if args.verbose:
			print("    There are", count, "connected components (CC)")
		CCs = (G.subgraph(CCs) for CCs in nx.connected_components(G)) # When using a 'networkx' version above 2.4
		#CCs = nx.connected_component_subgraphs(G)   # When using a 'networkx' version below 2.1
		c = 0
		with open(outNetworkCCs, "a") as outfile:
			for CC in CCs:
				c = c + 1
				if args.verbose:
					print("\r      Working on CC ", c, "/", count, sep="", end="")
				assort = nx.attribute_assortativity_coefficient(CC, "groups")
				if str(assort) == "nan":
					assort = 0
				attrib = nx.get_node_attributes(CC, "groups")
				attrib = attrib.values()
				attrib = set(attrib)
				if None in attrib:
					attrib = set(list(filter(None, attrib)))
				attrib = sorted(attrib)
				attrib = "|".join(attrib)			
				print((str(list(tmp.columns)[1]) + "\t" + str(c) + "\t" + str(assort) + "\t" + str(attrib)), file=outfile)
	if args.verbose:
		print("\n")
if args.verbose:
	print("Done")
