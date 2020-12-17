#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import re
import os

parser = argparse.ArgumentParser(description="Split a fasta file at given positions. Output files will be exported to the input file name followed by increasing integers.")

# Add the arguments to the parser
parser.add_argument("-f", "--file", dest="file_in", required=True,
                    help="A fasta file (.fasta, .fst or .fa)")

parser.add_argument("-p", "--positions", dest="position", required=True,
                    help="Positions to split the fasta file given in a string and separated by a '+' (i.e.; '1832+2204'). Each position should be the position at which each file will start. If you want to split from beginning to end of the fasta file, include a plus at the beginning and/or end of the string (i.e.; '+1832+2204+')")

parser.add_argument("-o", "--outputDir", dest="directory", required=False, default=None,
                    help="If selected, output files will be exported to the given directory.")

parser.add_argument("-r", "--remove", dest="remove", required=False, default=None, action="store_true",
                        help="If selected, will remove sequences composed only by gaps ('-').")

parser.add_argument("-u", "--unaligned", dest="unalign", required=False, default=None, action="store_true",
                        help="If selected, removes gaps ('-') in output files.")

args = parser.parse_args()

position = re.sub("^\+", "0+", args.position)
position = re.sub("\+$", "+end", position)

file_i = 0
positions = {}
for pos in position.split('+'):
    file_i += 1
    positions[pos] = file_i

if args.directory is not None:
    d = args.directory + "/"
    if not os.path.exists(args.directory):
        os.makedirs(args.directory)
else:
    d = ""

for i in list(range(0, len(positions)-1)):
    b=list(positions.keys())[i]
    e=list(positions.keys())[i+1]
    fileout = d + ".".join(args.file_in.split(".")[:-1]) + "_" + str(i+1) + "." + re.sub(".*\.", "", args.file_in)
    print("  Writing file ", int(i+1), " to '", fileout, "'", sep='')
    r = 0
    u = 0
    with open(fileout, "w") as outfile:
        for f in SeqIO.parse(open(args.file_in), "fasta"):
            if e == "end":
                 f.seq = f.seq[int(b):]
            else:
                f.seq = f.seq[int(b):(int(e))]
            if args.remove is not None:
                if list(set(f.seq)) != list('-'):
                    if args.unalign is not None:
                        f.seq = f.seq.ungap("-")
                        SeqIO.write(f, outfile, 'fasta')
                    else:
                        SeqIO.write(f, outfile, 'fasta')
                else:
                    r += 1
            else:
                if args.unalign is not None:
                    f.seq = f.seq.ungap("-")
                    SeqIO.write(f, outfile, 'fasta')
                else:
                    SeqIO.write(f, outfile, 'fasta')
            if args.remove is None:
                if list(set(f.seq)) == list(''):
                    u += 1
    if args.remove is not None:
        if r > 0: 
            print("   ", r, "sequence(s) contain only gaps ('-'), so they were removed.")
    if args.remove is None:
        if u > 0:
            print("   Warning! There are", u, "empty sequence(s). Consider using the '-r/--remove' option.")

