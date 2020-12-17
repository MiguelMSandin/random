#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import sys

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

out = {}
for i in list(range(0, len(files))):
    filen = list(files.keys())[i]
    tmp = files[filen]
    if i == 0:
        for key, value in tmp.items():
            out[key] = value
    else:
        for key, value in tmp.items():
            if key in out.keys():
                out[key] = (out[key] + value)
            elif key not in out.keys() and args.align is not None:
                r = len(out[list(out.keys())[1]])
                out[key] = ("-" * r + value)
            elif key not in out.keys():
                out[key] = value
        for key, value in out.items():
            if key not in tmp.keys():
                if args.align is not None:
                    r = len(tmp[list(tmp.keys())[1]])
                    out[key] = (value + "-" * r)

if args.drop is not None:
    d = len(out.keys())
    toRemove = set()
    for key in out.keys():
        for i in files:
            tmp = files[i]
            if key not in tmp.keys():
                toRemove.add(key)
    for i in toRemove:
        del out[i]

with open(args.file_out, "w") as outfile:
    for line in out:
        l = ">" + str(line) + "\n" + str(''.join(out[line])) + "\n"
        print(l, file=outfile)

print("  Final file contains", len(out), "sequences.")
if args.align is not None:
    tmp = out[list(out.keys())[1]]
    tmp = ''.join(tmp)
    tmp = len(tmp)    
    l = set()
    for s in out.values():
        l.add(len(s))
    if len(l) > 1:
        print("    Warning! sequences in input files were not aligned.")
    else:
        print("    And has", tmp, "aligned positions.")
if args.drop is not None:
    print("   ", d-len(out.keys()), "sequences were in the input files and are not in the final fasta.")

if args.drop is None and args.align is None:        
    l = set()
    for s in out.values():
        l.add(len(s))
    g = 0
    for val in out.values():
        if "-" in val:
            g += 1
    if g > 0 and len(l) > 0:
        print("  Warning!\n    You haven't selected any option to deal with gaps.\n    The final file has sequences of different length and with gaps.")
