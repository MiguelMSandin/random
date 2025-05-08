#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="From a network and a tab separated table with the nodes and attributes, summarises the network by counting connections between given attributes.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
					help="Network file. At least 2 columns for the given source and target nodes. The rest of the columns will be ignored.")

requiredArgs.add_argument("-t", "--table", dest="table", required=True,
					help="A tab delimited table with two columns: nodes and attributes. Other columns (if present) will be ignored.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
					help="The output file name. By default will replace the extension of the input net by '_summary.net'")

parser.add_argument("-n", "--headers", dest="headers", required=False, action="store_true",
					help="If selected, will assume the input network has headers, and ignore the first entry.")

parser.add_argument("-a", "--addHeaders", dest="addHeaders", required=False, action="store_true",
					help="If selected, will add the following headers to the net: 'source target count'.")

parser.add_argument("-s", "--standardization", dest="standardization", required=False, default=None,
					help="Whether the output count is standardize or not. Available standardizations are 'total', 'max', 'log' and 'hellinger'. For details on the standardization see the R function 'decostand' from the package 'vegan'. By default, no standardization is applied.")

parser.add_argument("-b", "--base", dest="base", required=False, default=None,
					help="The base of the logarithmic standardization when 'log' is selected. By default exp(1)")

args = parser.parse_args()

if args.file_out is None:
	import re
	outfile = re.sub("\\.[^\\.]+$", "_summary.net", args.file_in)
else:
	outfile = args.file_out

if args.standardization is not None:
	if args.standardization not in ['total', 'max', 'log', 'hellinger']:
		print("  Error! Transformation has to be one of 'total', 'max', 'log', 'hellinger' or 'None'.")
		print("Exiting")
		import sys
		sys.exit(0)
	import math

# Reading table
print("  Reading table")
table = {}
for line in open(args.table):
	line = line.strip().split()
	table[line[0]] = line[1]

# Reading table
print("  Reading and summarising network")
i = 0
out = {}
for line in open(args.file_in):
	i += 1
	if args.headers and i == 1:
		print("    Ignoring network headers")
	else:
		line = line.strip().split()
		a1 = table[line[0]]
		a2 = table[line[1]]
		tmp = [a1, a2]
		tmp.sort()
		attr = str(tmp[0]) + "\t" + str(tmp[1])
		if attr in out.keys():
			out[attr] += 1
		else:
			out[attr] = 1

# If selected, set parameters for standardization
if args.standardization is not None:
	print("  Applying a '", args.standardization, "' standardization", sep="")
	if args.standardization == 'total' or args.standardization == 'hellinger':
		total = 0
		for val in out.values():
			total += val
	if args.standardization == 'max':
		maxv = 0
		for val in out.values():
			if val > maxv:
				maxv = val
	if args.standardization == 'log':
		if args.base is None:
			base = math.exp(1)
		else:
			base = args.base

# Exporting summarised network
print("  Exporting summarised table to:", outfile)
print("    Edges in network:\t", i, sep="")
print("    Edges summarised:\t", len(out), sep="")
with open(outfile, "w") as outFile:
	if args.addHeaders:
		print("source\ttarget\tcount", file=outFile)
	for key, val in out.items():
		if args.standardization == 'total':
			valo = val / total
		elif args.standardization == 'max':
			valo = val / maxv
		elif args.standardization == 'log':
			valo = math.log(val, base)
		elif args.standardization == 'hellinger':
			valo = math.sqrt(val / total)
		else:
			valo = val
		print(key + "\t" + str(valo), file=outFile)

print("Done")
