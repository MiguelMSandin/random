#!/usr/bin/env python3

import argparse
from collections import defaultdict
import re

parser = argparse.ArgumentParser(description="Creates a network file from a similarity blast output.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
					help="Input file. This assumes a file with the following columns: 'qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen'.")

requiredArgs.add_argument("-i", "--identity", dest="identity", required=True, nargs='+',
					help="Identity threshold(s) to establish a connection between nodes.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
					help="Output file. Returns the network file. Files will be exported to the output file name ('output') followed by the selected thresholds and the '.net' extension (i.e.; 'output_80.net'). By default will add '_ID.net' to the input file.")

parser.add_argument("-c", "--cover", dest="cover", required=False, default=None,
					help="Applies a minimal cover threshold to establish a connection between nodes, in addition to the identity and evalue, if selected, thresholds.")

parser.add_argument("-e", "--evalue", dest="evalue", required=False, default=None,
					help="Applies a maximal evalue threshold to establish a connection between nodes, in addition to the identity and cover, if selected, thresholds.")

parser.add_argument("-I", "--identityPosition", dest="identityPosition", required=False, default=4, type=int,
					help="The positional column of the identity value. Default=4 (in the fourth column).")

parser.add_argument("-k", "--clean", dest="cleaning", required=False, action="store_true", default=None,
					help="If selected, the input file is cleaned from reciprocal (A-B = B-A) and/or equal (A-A) hits before creating the network file.")

parser.add_argument("-a", "--addHeaders", dest="addHeaders", required=False, action="store_true",
					help="If selected, will add the following headers to the net: 'source target id'.")

args = parser.parse_args()

if args.file_out is None:
	out = re.sub("\\.[^\\.]+$", "", args.file_in)
else:
	out = args.file_out

# Clean the file if needed
if args.cleaning is not None:
	print("  Cleaning the file:")
	print("  Setting hit names")
	clean = {}
	removed = 0
	accepted = 0
	selfHit = 0
	for line in open(args.file_in):
		line = line.strip().split()
		seq1 = line[0]
		seq2 = line[1]
		if seq1 == seq2:
			selfHit = selfHit + 1
		hit = str(seq1) + "\t" + str(seq2)
		hitr =str(seq2) + "\t" + str(seq1)
		if hit in clean:
			removed = removed + 1
		elif hitr in clean:
			removed = removed + 1
		else:
			values = str(line[2]) + "\t" + str(line[3]) + "\t" + str(line[4]) + "\t" + str(line[5]) + "\t" + str(line[6]) + "\t" + str(line[7]) + "\t" + str(line[8]) + "\t" + str(line[9]) + "\t" + str(line[10])
			clean[hit] = values
			accepted = accepted + 1
	print("  Writing cleaned file with '", accepted, "' hits", sep="")
	cleanFile = re.sub(r"\..*$", "_clean.similarities", args.file_in)
	with open(cleanFile, "w") as outfile:
		for hit in list(clean.keys()):
			print(hit + "\t" + clean[hit], file=outfile)
	print("  Removed hits: ", removed)
	print("  Self matching:", selfHit)
	print("  Cleaned file exported to: '", re.sub(r"\..*$", "_clean.similarities", args.file_in), "'", sep="")

# Setting name of input file
if args.cleaning:
	input_file = re.sub(r"\..*$", "_clean.similarities", args.file_in)
else:
	input_file = args.file_in

# Now looping through the different identity thresholds (if more than one) and selecting nodes
for i in args.identity:
	if args.cover is not None:
		c=", '" + str(args.cover) + "' minimum cover"
		C="_" + "c" + str(args.cover)
	else:
		c=""
		C=""
	if args.evalue is not None:
		e=", '" + str(args.evalue) + "' maximum evalue"
		E="_" + "e" + str(args.args)
	else:
		e=""
		E=""
	print("  Creating a network file with '", i, "%' similarity threshold", c, e, sep="")
	tmp = out + "_" + "i" + i + C + E + ".net"
	countin = 0
	countout = 0
	with open(tmp, "w") as outfile:
		if args.addHeaders:
			print("source\ttarget\tid", file=outfile)
		for line in open(input_file):
			countin += 1
			data = line[:-1].split("\t")
			seq1 = data[0]
			seq2 = data[1]
			seqid = data[args.identityPosition-1]
			if float(seqid) >= float(i):
				statement1=True
			else:
				statement1=False
			if args.cover is not None:
				cover1 = 100.0 * (int(data[6])-int(data[5])) / int(data[7])
				cover2 = 100.0 * (int(data[9])-int(data[8])) / int(data[10])
				if (cover1 >= int(args.cover)) and (cover2 >= float(args.cover)):
					statement2=True
				else:
					statement2=False
			else:
				statement2=None
			if args.evalue is not None:
				ev = data[2]
				if (float(ev) <= float(args.evalue)):
					statement3=True
				else:
					statement3=False
			else:
				statement3=None
			if statement1 and (statement2 or statement2==None) and (statement3 or statement3==None):
				countout += 1
				lineOut = str(seq1) + "\t" + str(seq2) + "\t" + str(seqid)
				print(lineOut, file=outfile)
	print("  Network file exported to: '", tmp, "'", sep="")
	print("    Nodes in:\t", countin, sep="")
	print("    Nodes out:\t", countout, "\t(", round(countout/countin*100, 1), "%)", sep="")

print("Done")
