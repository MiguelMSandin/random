#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import re
import sys


parser = argparse.ArgumentParser(description="Remove sequences in a fasta file that are duplicated and (optional) that are shorter than 'length' bp.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input fasta file.")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
                    help="Output file. By default will add '_clean' to file name before the extension.")

parser.add_argument("-l", "--length", dest="length", required=False, default=None, type=float, action="store",
                    help="If selected, will remove sequences shorter than 'length' bp.")

parser.add_argument("-d", "--duplicates", dest="duplicates", required=False, default=None, action="store_true",
                    help="If selected, will remove duplicated sequences.")

parser.add_argument("-i", "--headers", dest="headers", required=False, default=None, action="store_true",
                    help="If selected, will only check for duplicate sequence names.")

parser.add_argument("-s", "--sequences", dest="sequences", required=False, default=None, action="store_true",
                    help="If selected, will only check for duplicate sequences.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, default=None, action="store_false",
                    help="If selected, will NOT print information to the console.")

args = parser.parse_args()

# Set output name
if args.file_out is None:
	file_out = re.sub("\\.[^\\.]+$", "_clean", args.file_in) + re.sub(".*\.", ".", args.file_in)
else:
	file_out = args.file_out

# Decide if duplicates are going to be removed
if args.sequences is not None and args.headers is not None:
	duplicates = True
elif args.duplicates is not None:
	duplicates = True
else:
	duplicates = False

# Start cleaning
seq_in = 0
seq_out = 0
with open(file_out, "w") as outfile:
	seqsid = set()
	seqs = set()
	for i in SeqIO.parse(open(args.file_in), "fasta"):
		seq_in += 1
		seqid_i=i.id
		seq_i=str(i.seq)
		writing = False
		if args.headers is not None and args.sequences is None:
			if seqid_i not in seqsid:
				writing = True
		elif args.sequences is not None and args.headers is None:
			if seq_i not in seqs:
				writing = True
		if duplicates:
			if (seqid_i not in seqsid) and (seq_i not in seqs):
				writing = True
		if args.length is not None:
			tmp = re.sub("-", "", seq_i)
			if len(tmp) >= args.length:
				writing = True
		if writing:
			SeqIO.write([i], outfile, "fasta")
			seq_out += 1
		seqsid.add(seqid_i)
		seqs.add(seq_i)

if not args.verbose:
	print("  Sequences in: ", seq_in)
	print("  Sequences out:", seq_out)
 
