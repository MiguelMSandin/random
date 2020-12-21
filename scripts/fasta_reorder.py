#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Phylo

parser = argparse.ArgumentParser(description="Reorders a fasta file based on a rooted tree.")

# Add the arguments to the parser
parser.add_argument("-f", "--fasta", dest="file_in", required=True,
                    help="Fasta files.")

parser.add_argument("-t", "--tree", dest="tree", required=True,
                        help="Rooted tree file.")

parser.add_argument("-F", "--format", dest="format", required=False, default="nexus",
                        help="Tree file format, defaul='nexus'")

parser.add_argument("-o", "--output", dest="file_out", required=True,
                    help="The output fasta file name.")

parser.add_argument("-d", "--drawTree", dest="draw", required=False, action="store_true", default=None,
                    help="If selected, draws the tree in the console while reordering the fasta file.")

args = parser.parse_args()

print("  Reading tree")
T = Phylo.read(args.tree, args.format)
if args.draw is not None:
    Phylo.draw_ascii(T)

print("  Reading fasta")
fasta = {}
for line in SeqIO.parse(open(args.file_in), "fasta"):
    fasta[line.id] = str(line.seq)

print("  Reordering and writing fasta")
with open(args.file_out, "w") as outfile:
    for line in  T.get_terminals():
        sp = line.name
        sp = sp.strip("'")
        print(">" + str(sp) + "\n" + str(fasta[sp]), file=outfile)

print("Done")
