#!/usr/bin/env python3

import argparse
import re

parser = argparse.ArgumentParser(description="Changes the names of a fasta file given a tab delimited table with the old names in one column and the new names in the second column.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input fasta file.")

requiredArgs.add_argument("-l", "--list", dest="table", required=True,
                    help="A tab delimited table with the old names in one column and the new names in the second column.")

parser.add_argument("-o", "--output", dest="file_out", required=False,
                    help="Output fasta file. By default will add '_newName' to the input file name before the extension.")

parser.add_argument("-r", "--order", dest="order", required=False, action="store_false",
					help="If selected, will change the names according to the order of the given list. In this case, the list must have only one column corresponding to the new names.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Setting variables --------------------------------------------------------------------------------
if args.file_out is None:
	outFile = re.sub("\\.[^\\.]+$", "_newName.", args.file_in) + re.sub(".*\\.", "", args.file_in)
else:
	outFile = args.file_out

# Reading table ------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading table:", args.table)
if args.order:
	names = {}
	for line in open(args.table):
		tmp = line.strip().split()
		oldname = tmp[0]
		newname = tmp[1]
		newname = re.sub("\n", "", newname)
		names[oldname] = newname
else:
	names = list()
	countl = 0
	for line in open(args.table):
		countl += 1
		name = line.strip().split()
		name = name[0]
		name = re.sub("\n", "", name)
		names.append(name)

# Reading table ------------------------------------------------------------------------------------
if args.order:
	if args.verbose:
		print("  Changing names")
	notFound = list()
	with open(outFile, "w") as outfile:
		for line in open(args.file_in):
			if ">" in line:
				oldname = re.sub(r"\>", "", line)
				oldname = re.sub("\n", "", oldname)
				if oldname in names:
					newname = names[oldname]
					print(">" + newname, end="\n", file=outfile)
				else:
					notFound.append(oldname)
					print(">" + oldname, end="\n", file=outfile)
			else:
				print(line, end="", file=outfile)
	if len(notFound) > 0:
		print("  The following tips were not found in the list, and therefore not renamed:")
		for name in notFound:
			print("   ", name, end="")
else:
	if args.verbose:
		print("  Changing names in the given order")
	with open(outFile, "w") as outfile:
		countf = 0
		for line in open(args.file_in):
			if ">" in line:
				print(">" + names[countf], end="\n", file=outfile)
				countf += 1
			else:
				print(line, end="", file=outfile)
	if countl != countf:
		print("  Warning! Different number of sequences (", countf, ") and names (", countl, ") provided.", sep="")

if args.verbose:
	print("Done")
