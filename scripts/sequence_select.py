#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import re
import sys

parser = argparse.ArgumentParser(description="Select sequence in a fasta file either from a list or that match a pattern and extract them or remove them.")

parser.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input fasta file.")

parser.add_argument("-o", "--output", dest="file_out", required=True,
                    help="Output file.")

parser.add_argument("-l", "--list", dest="listSeq", required=False, default=None,
                    help="List of sequences to be selected. This must be a different file with each sequence name in a different line.")

parser.add_argument("-p", "--pattern", dest="pattern", required=False, default=None,
                    help="Pattern to be matched for selection of the sequences.")

parser.add_argument("-a", "--action", dest="action", required=True, default='keep', choices=['k', 'keep', 'r', 'remove'],
                    help="Either 'remove' or 'keep' (default), if you want to remove or export (keep) selected sequences. Their initials are also working.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, default=None, action="store_true",
                    help="If selected, will print information of sequences selected.")

args = parser.parse_args()
   
if args.listSeq is None and args.pattern is None:
    print("Error: You should specify either a list (-l/--list) or a pattern (-p/--pattern) so I can select sequences...")
    sys.exit(1)

if args.listSeq is not None:
    listSeq = [line.strip() for line in open(args.listSeq)]

seq_in = 0
seq_out = 0
with open(args.file_out, "w") as outfile:
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
                if args.pattern not in i.id:
                    SeqIO.write([i], outfile, "fasta")
                    seq_out += 1
        if args.action == "keep" or args.action == "k":
            if args.listSeq is not None:
                if i.id in listSeq:
                    SeqIO.write([i], outfile, "fasta")
                    seq_out += 1
            if args.pattern is not None:
                if args.pattern in i.id:
                    SeqIO.write([i], outfile, "fasta")
                    seq_out += 1

if args.verbose is not None:
    print("")
    if args.listSeq is not None:
        missing = set()
        for l in listSeq:
            if l not in seqsid:
                missing.add(l)
        if len(missing) > 0:
            print("  Sequences present in the list and not found in the input fasta file:\n", "\n ".join(missing))
    print("  Sequences in: ", seq_in)
    print("  Sequences out:", seq_out)
        
