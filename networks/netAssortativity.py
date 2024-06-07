#!/usr/bin/env python3

import argparse
import networkx as nx
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Analyzes the assortativity of groups from a network file(s).")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True, nargs="+",
					help="Input file. A network file with three columns: 'qseqid sseqid id'.")

requiredArgs.add_argument("-a", "--attributes", dest="file_attr", nargs="+", required=True,
					help="Attributes file(s). A tab separated file(s) with the nodes in the first column and the rest of the groups in the following columns matching the given node. The order of the files should match that given for the networks. Note: Remember to add the column names, so it takes that name to sort out the output.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
					help="Output file. Returns a tab delimited file with the following output: the name of the file, the connected component, the number of nodes, the number of edges, the assortativity, and the unique name of states for the give connected component. By default, will print the information to the console. If the provided file exists, it will append the information to it.")

parser.add_argument("-i", "--ignore", dest="ignore", required=False, default=None, type=float,
					help="If selected, will ignore connected components with less nodes than the given value. If (0-1] will take the percentage of the whole graph.")

parser.add_argument("-c", "--headers", dest="headers", required=False, action="store_false",
					help="If selected, will skip the first row of the network interpreted as headers.")

parser.add_argument("-n", "--names", dest="names", required=False, action="store_true",
					help="If selected, will print to the given output name the column names. If output name not selected, this command will be ignored.")

parser.add_argument("-d", "--overwrite", dest="delete", required=False, action="store_true",
					help="If selected, will overwrite the output file if exists.")

parser.add_argument("-s", "--simple", dest="simple", required=False, action="store_true",
					help="If selected, will just print the assortativity value to the console")

parser.add_argument("-r", "--randomize", dest="random", required=False, default=None, type=int,
					help="If selected, will randomize the categories of the given attribute N times to test for deviation of random assortativities and estimate the probabilities of debiating from random assortativity with a one-sample T-test.")

args = parser.parse_args()

if len(args.file_in) != len(args.file_attr):
	print("  Error: The number of networks (", len(args.file_in), ") is different than the number of atrtibute files (", len(args.file_attr), ") provided\nExiting", sep="")
	import sys
	sys.exit(1)

# Define a function to check if a given net is a clique
def is_clique(G):
    n = len(list(G.nodes()))
    e = len(list(G.edges()))
    return e == n*(n-1)/2

files = {}
for i in range(len(args.file_in)):
	net = args.file_in[i]
	atr = args.file_attr[i]
	files[net] = atr

if args.ignore is not None:
	if args.ignore < 0:
		print("  Warning! CCs smaller than 0 can't exist. Parameter 'ignore' set to 0, which won't have any effect")
		ignore = 0
	else:
		ignore = args.ignore
else:
	ignore = 0

if args.random is not None:
	import random
	import statistics
	from scipy import stats

if (args.file_out is None) and (args.simple == False):
	if args.random is not None:
		print("file_net\tfile_attribute\tCC\tnodes\tedges\tassortativity\trandom\tsignificance\tattribute\tstates")
	else:
		print("file_net\tfile_attribute\tCC\tnodes\tedges\tassortativity\tattribute\tstates")
elif args.file_out is not None:
	out = args.file_out
	if args.delete:
		import os.path
		if os.path.exists(out):
			os.remove(out)
	if args.names:
		with open(out, "a") as outfile:
			if args.random is not None:
				print("file_net\tfile_attribute\tCC\tnodes\tedges\tassortativity\trandom\tsignificance\tattribute\tstates", file=outfile)
			else:
				print("file_net\tfile_attribute\tCC\tnodes\tedges\tassortativity\tattribute\tstates", file=outfile)

for net, atr in files.items():
	if (args.file_out is None) and (args.simple == False):
		print(str(net), end="")
	# Reading network
	if args.headers:
		G=nx.read_edgelist(net, delimiter="\t", data=(("id",float),))
	else:
		with open(net, 'rb') as neti:
			next(neti, '')
			G = nx.read_edgelist(neti, delimiter='\t', data=(("id",float),))
		
	# Reading attributes
	attributesTable = pd.read_csv(atr, sep="\t")
	
	# Now loop through the attributes file to work on the different attributes
	for i in list(range(1, len(attributesTable.columns))):
		column = pd.concat([attributesTable.iloc[:,0], attributesTable.iloc[:,i]], axis=1)
		attr = {}
		for line in column.index:
			key = column.iloc[:,0][line]
			val = column.iloc[:,1][line]
			if str(val) == "nan":
				val = "NA"
			attr[key] = val
		
		# Count number of nodes and edges
		nodes=list(G.nodes())
		edges=list(G.edges())
		
		# Now matching order of every node with the attribute
		group = {}
		attributes = []
		for node in nodes:
			tmp = attr.get(node)
			group[node] = tmp
			if tmp not in attributes:
				if tmp is not None:
					attributes.append(tmp)
				if (tmp is None) & ("NA" not in attributes):
					attributes.append("NA")
		attributes = sorted(attributes)
		
		# Calculate assortativity
		if (len(attributes) == 1) or (is_clique(G)):
			assort = "nan"
		else:
			# Adding the attributes to the network
			nx.set_node_attributes(G, group, "groups")
			assort = nx.attribute_assortativity_coefficient(G, "groups")
			assort = round(assort, 12)
		
		if (args.random is not None) and (args.simple == False):
			if (len(attributes) == 1) or (is_clique(G)):
				ttest = "nan"
				ranMean = "nan"
			else:
				assortr = []
				for i in range(args.random):
					group = {}
					for node in nodes:
						tmp = random.sample(attributes, k=1)[0]
						group[node] = tmp
					nx.set_node_attributes(G, group, "groups")
					assorti = nx.attribute_assortativity_coefficient(G, "groups")
					assortr.append(assorti)
				stats.ttest_1samp(assortr, assort)
				ttest = stats.ttest_1samp(assortr, assort)[1]
				ttest = round(ttest, 12)
				ranMean = statistics.mean(assortr)
				ranMean = round(ranMean, 12)
			
		# And print information
		attributes = "|".join(attributes)
		attribute = str(list(column.columns)[1])
		if (args.file_out is None) and (args.simple == False):
			if args.random is not None:
				print("\r", str(net) + "\t" + str(atr) + "\tGraph\t" + str(len(nodes)) + "\t" + str(len(edges)) + "\t" + str(assort) + "\t" + str(ranMean) + "\t" + str(ttest) + "\t" + str(attribute) + "\t" + str(attributes), sep="")
			else:
				print("\r", str(net) + "\t" + str(atr) + "\tGraph\t" + str(len(nodes)) + "\t" + str(len(edges)) + "\t" + str(assort) + "\t" + str(attribute) + "\t" + str(attributes), sep="")
		elif args.simple:
			print(str(assort))
		else:
			with open(out, "a") as outfile:
				if args.random is not None:
					print(str(net) + "\t" + str(atr) + "\tGraph\t" + str(len(nodes)) + "\t" + str(len(edges)) + "\t" + str(assort) + "\t" + str(ranMean) + "\t" + str(ttest) + "\t" + str(attribute) + "\t" + str(attributes), file=outfile)
				else:
					print(str(net) + "\t" + str(atr) + "\tGraph\t" + str(len(nodes)) + "\t" + str(len(edges)) + "\t" + str(assort) + "\t" + str(attribute) + "\t" + str(attributes), file=outfile)
		
		# Check if there are more than one connected commponents (CCs)
		if args.simple == False:
			CCs = (G.subgraph(CCs) for CCs in nx.connected_components(G))
			count = 0
			for CC in CCs:
				count += 1
			# if there is more than one, repeat previous steps for each CC
			if count > 1:
				j = 0
				CCs = (G.subgraph(CCs) for CCs in nx.connected_components(G))
				for CC in CCs:
					j += 1
					nodes=list(CC.nodes())
					n = len(nodes)
					edges=list(CC.edges())
					e = len(edges)
					if (ignore > 0) & (ignore <= 1):
						ignore = n * ignore
					if n >= ignore:
						# Now matching order of every node with the attribute
						group = {}
						attributes = []
						for node in nodes:
							tmp = attr.get(node)
							group[node] = tmp
							if tmp not in attributes:
								if tmp is not None:
									attributes.append(tmp)
								if (tmp is None) & ("NA" not in attributes):
									attributes.append("NA")
						attributes = sorted(attributes)
						# Calculate assortativity
						if (len(attributes) == 1) or (is_clique(CC)):
							assort = "nan"
						else:
							nx.set_node_attributes(CC, group, "groups")
							assort = nx.attribute_assortativity_coefficient(CC, "groups")
							assort = round(assort, 12)
						if args.random is not None:
							if (len(attributes) == 1) or (is_clique(CC)):
								ranMean = "nan"
								ttest = "nan"
							else:
								assortr = []
								for i in range(100):
									group = {}
									for node in nodes:
										tmp = random.sample(attributes, k=1)[0]
										group[node] = tmp
									nx.set_node_attributes(CC, group, "groups")
									assorti = nx.attribute_assortativity_coefficient(CC, "groups")
									assortr.append(assorti)
								stats.ttest_1samp(assortr, assort)
								ttest = stats.ttest_1samp(assortr, assort)[1]
								ttest = round(ttest, 12)
								ranMean = statistics.mean(assortr)
								ranMean = round(ranMean, 12)
						# And print the information
						attributes = "|".join(attributes)
						if args.file_out is None:
							if args.random is not None:
								print("\r", str(net) + "\t" + str(atr) + "\t" + str(j) + "\t" + str(len(nodes)) + "\t" + str(len(edges)) + "\t" + str(assort) + "\t" + str(ranMean) + "\t" + str(ttest) + "\t" + str(attribute) + "\t" + str(attributes), sep="")
							else:
								print("\r", str(net) + "\t" + str(atr) + "\t" + str(j) + "\t" + str(len(nodes)) + "\t" + str(len(edges)) + "\t" + str(assort) + "\t" + str(attribute) + "\t" + str(attributes), sep="")
						else:
							with open(out, "a") as outfile:
								if args.random is not None:
									print(str(net) + "\t" + str(atr) + "\t" + str(j) + "\t" + str(len(nodes)) + "\t" + str(len(edges)) + "\t" + str(assort) + "\t" + str(ranMean) + "\t" + str(ttest) + "\t" + str(attribute) + "\t" + str(attributes), file=outfile)
								else:
									print(str(net) + "\t" + str(atr) + "\t" + str(j) + "\t" + str(len(nodes)) + "\t" + str(len(edges)) + "\t" + str(assort) + "\t" + str(attribute) + "\t" + str(attributes), file=outfile)
