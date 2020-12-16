#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import sys
#import os

parser = argparse.ArgumentParser(description="Concatenate multiple fasta files. Bear in mind that the sequence names should be exactly the same in every file.")

# Add the arguments to the parser
parser.add_argument("-f", "--files", dest="files_in", nargs='+', required=True,
                    help="Fasta files to be concatenated in the given order.")

parser.add_argument("-o", "--output", dest="file_out", required=True,
                    help="The output file name.")

parser.add_argument("-a", "--align", dest="align", required=False, default=None, action="store_true",
                        help="If selected, will replace by gaps ('-') those sequences not present in a file. This understands you are providing aligned fasta files.")

parser.add_argument("-d", "--drop", dest="drop", required=False, default=None, action="store_true",
                        help="If selected, will only print sequences present in all files.")

args = parser.parse_args()

if args.align is not None and args.drop is not None:
    print("\nError: Either you want to keep the sequences or not, but '-a/--align' and '-d/--drop' arguments are doing kind of opposite stuff.\n")
    sys.exit(1)

files = {}
for filei in args.files_in:
    f = {}
    for line in SeqIO.parse(open(filei), "fasta"):
        f[line.id] = str(line.seq)
    files[filei] = f

out = files[list(files.keys())[0]]

d = 0
for i in list(range(1, len(files))):
    filen = list(files.keys())[i]
    tmp = files[filen]
    for key, value in tmp.items():
        if key in out:
            out[key] = [out[key], value]
        elif key in out and args.drop is not None:
            d += 1
            del out[key]
        elif key not in out and args.align is not None:
            r = len(out[list(out.keys())[1]])
            out[key] = ["-" * r, value]
    for key, value in out.items():
        if key not in tmp.keys():
            r = len(tmp[list(tmp.keys())[1]])
            out[key] = [value, "-" * r]

with open(args.file_out, "w") as outfile:
    for line in out:
        l = ">" + str(line) + "\n" + str(''.join(out[line])) + "\n"
        print(l, file=outfile)

print("  Final file contains", len(out), "sequences.")
if args.align is not None:
    tmp = out[list(out.keys())[1]]
    tmp = ''.join(tmp)
    tmp = len(tmp)
    print("    And has", tmp, "aligned positions.")
if args.drop is not None:
    print("   ", d, "sequences were in the input files and are not in the final fasta.")
    
    
