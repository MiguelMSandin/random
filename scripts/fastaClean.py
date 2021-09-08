#!/usr/bin/env python3

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Remove sequences in a fasta file that are duplicated and (optional) that are shorter than 'length' bp.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input fasta file.")

requiredArgs.add_argument("-o", "--output", dest="file_out", required=True,
                    help="Output file.")

parser.add_argument("-i", "--headers", dest="headers", required=False, default=None, action="store_true",
                    help="If selected, will only check for duplicate sequence names.")

parser.add_argument("-s", "--sequences", dest="sequences", required=False, default=None, action="store_true",
                    help="If selected, will only check for duplicate sequences.")

parser.add_argument("-l", "--length", dest="length", required=False, default=None, type=float, action="store",
                    help="If selected, will remove sequences shorter than 'length' bp.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, default=None, action="store_true",
                    help="If selected, will print information to the console.")

args = parser.parse_args()

seq_in = 0
seq_out = 0
with open(args.file_out, "w") as outfile:
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
		else:
			if (seqid_i not in seqsid) and (seq_i not in seqs):
				writing = True
		if args.length is not None:
			if len(seq_i) < args.length:
				writing = False
		if writing:
			SeqIO.write([i], outfile, "fasta")
			seq_out += 1
		seqsid.add(seqid_i)
		seqs.add(seq_i)

if args.verbose is not None:
	print("  Sequences in: ", seq_in)
	print("  Sequences out:", seq_out)
