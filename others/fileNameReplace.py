#!/usr/bin/env python3

import argparse
import re

parser = argparse.ArgumentParser(description="From a text file a tab separated table, will search for names in the table to replace them.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="input", required=True,
						  help="A text file.")

requiredArgs.add_argument("-l", "--list", dest="table", required=True,
					help="A tab delimited file with the names to be renamed in one column and the new name in a second column.")

parser.add_argument("-o", "--output", dest="output", required=False, default=None,
					help="The output file name. By default will add to the file name '_renamed' and the original extension.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Setting variables --------------------------------------------------------------------------------

# Output file
if args.output is None:
	outFile = re.sub("\\.[^\\.]+$", "_renamed.", args.input) + re.sub(".*\\.", "", args.input)
else:
	outFile = args.output

# Reading file -------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading list file:", args.table)
names = {}
for line in open(args.table):
	tmp = line.strip().split()
	names[tmp[0]]=tmp[1]

with open(outFile, "w") as outfile:
	if args.verbose:
		print("  Replacing text file:", args.input)
	for line in open(args.input):
		line = line.strip().split()
		out = list()
		for i in line:
			if i in names:
				out.append(names[i])
			else:
				out.append(i)
		out = " ".join(out)
		print(out, file=outfile)

if args.verbose:
	print("  New text file written to:" outFile)
	print("Done")
