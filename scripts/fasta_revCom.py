#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import sys

parser = argparse.ArgumentParser(description="Reverse and(/or) complement sequeces in a fasta file.")

parser.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input fasta file.")

parser.add_argument("-o", "--output", dest="file_out", required=True,
                    help="Output fasta file.")

parser.add_argument("-c", "--complement", dest="complement", required=False, default=None, action="store_true",
                    help="If selected, returns only the complement sequences.")

parser.add_argument("-r", "--reverse", dest="reverse", required=False, default=None, action="store_true",
                    help="If selected, returns only the reverse sequences.")

args = parser.parse_args()

if args.complement is not None and args.reverse is not None:
    print("\nError: Either you want the reverse or the complement sequences.\n  Please select either '-c/--complement', '-r/--reverse' or none, but not the two options.\n")
    sys.exit(1)

with open(args.file_out, "w") as outfile:
    for i in SeqIO.parse(open(args.file_in), "fasta"):
        sequence = i.seq
        name = i.id
        if args.complement is None and args.reverse is None:
            sequence = sequence.reverse_complement()
        if args.complement is not None:
            sequence = sequence.complement()
        if args.reverse is not None:
            sequence = sequence[::-1]
        print(">" + str(name) + "\n" + str(sequence) + "\n", file=outfile)
