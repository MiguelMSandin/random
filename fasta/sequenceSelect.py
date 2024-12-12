#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import re

parser = argparse.ArgumentParser(description="Select sequence in a fasta file either from a list or that match a pattern and extract them or remove them.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
						  help="Input fasta file. Avoid spaces in the sequence names.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
					help="Output file. By default will add '_selected' to the input file name. If the file already exists, sequences will be appended at the end of the file.")

parser.add_argument("-l", "--list", dest="listSeq", required=False, default=None,
					help="List of sequences to be selected. This must be a different file with each sequence name in a different line.")

parser.add_argument("-p", "--pattern", dest="pattern", required=False, default=None, nargs="+",
					help="Pattern(s) to be matched for selection of the sequences. When using this option in combination with '--action keep', it might be faster to use 'grep -A 1 PATTERN FILE_IN > FILE_OUT' if each sequence from the input file is in a single line.")

parser.add_argument("-a", "--action", dest="action", required=False, default='keep',
					choices=['k', 'keep', 'r', 'remove'],
					help="Either 'remove' or 'keep' (default), if you want to remove ('remove') or export ('keep') selected sequences. Their initials are also working.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

if args.listSeq is None and args.pattern is None:
	import sys
	print("Error: You should specify either a list (-l/--list) or a pattern (-p/--pattern) so I can select sequences...")
	sys.exit(1)

# Output file
if args.file_out is None:
	outFile = re.sub("\\.[^\\.]+$", "_selected", args.file_in) + re.sub(".*\\.", ".", args.file_in)
else:
	outFile = args.file_out

if args.listSeq is not None:
	listSeq = [line.strip() for line in open(args.listSeq)]

seq_in = 0
seq_out = 0
with open(outFile, "a") as outfile:
	seqsid = set()
	for i in SeqIO.parse(open(args.file_in), "fasta"):
		seq_in += 1
		seqsid.add(i.id)
		if args.action == "remove" or args.action == "r":
			if args.listSeq is not None:
				if i.id not in listSeq:
					SeqIO.write([i], outfile, "fasta")
					seq_out += 1
			if args.pattern is not None:
				for pattern in args.pattern:
					if pattern not in i.id:
						SeqIO.write([i], outfile, "fasta")
						seq_out += 1
		if args.action == "keep" or args.action == "k":
			if args.listSeq is not None:
				if i.id in listSeq:
					SeqIO.write([i], outfile, "fasta")
					seq_out += 1
			if args.pattern is not None:
				for pattern in args.pattern:
					if pattern in i.id:
						SeqIO.write([i], outfile, "fasta")
						seq_out += 1

if args.verbose:
	if args.listSeq is not None:
		missing = set()
		for l in listSeq:
			if l not in seqsid:
				missing.add(l)
		if len(missing) > 0:
			print("  The following sequences from the list are not found in the input fasta file:\n", "\n ".join(missing))
	print("  Sequences in: ", seq_in)
	print("  Sequences out:", seq_out)
