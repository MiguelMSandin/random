#!/usr/bin/env python3

import argparse
from Bio import AlignIO
import re
import sys
import math

parser = argparse.ArgumentParser(description="Calculates Shannon entropy, richness, unique bases, number of repetitions, the alignment coverage and/or the running mean of the Shannon entropy (mean shannon entropy at given window) at every position in an aligned fasta file.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="fastaFile", required=True,
                    help="An aligned fasta file.")

parser.add_argument("-o", "--output", dest="outFile", required=False, default=None,
                    help="The output name of the alignment information file. By default will add '_positionInfo.tsv' to file name excluding the extension.")

parser.add_argument("-c", "--columns", dest="fields", required=False, default=None,
                    help="Values or fields to be computed. By default: 'position+shannon+richness+unique+repetitions+cover+smooth'.")

parser.add_argument("-b", "--base", dest="base", required=False, action='store', type=int, default=math.e,
                    help="Base of the Shannon entropy. Default: 'e'.")

parser.add_argument("-r", "--range", dest="rangeSmooth", required=False, action='store', type=int, default=50,
                    help="The smooth value. Calculates the average Shannon entropy values at each position in a given window (position +- range; window=2*range). Default: 50.")

parser.add_argument("-g", "--gaps", dest="clean", required=False, default=None, action="store_true",
                        help="If selected, will also compute values of Shannon entropy removing gaps (-).")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, default=None, action="store_true",
                    help="If selected, will print information in the console.")

args = parser.parse_args()

# __________________________________________________________________________________________________
if args.verbose:
	print("  Reading fasta...", end="")
fasta =  AlignIO.read(args.fastaFile, "fasta")

length = set()
for seq in fasta:
	length.add(len(seq))

seqs = len(fasta)

if len(length) != 1:
	print("\nError: Input file is not aligned.\nExiting\n")
	sys.exit(1)
else:
	length = int(list(length)[0])

# __________________________________________________________________________________________________
if args.verbose:
	print("\r  Reading fasta and asigning positions...")
nucleotides = {}
for position in range(length):
	p = position +1
	if args.verbose:
		print("\r    ", p, "/", length, sep="", end="")
	nucleotides[p] = list(fasta[:, position])

if args.verbose:
	print("\n    Fasta file contains '", seqs, "' sequences and '", length, "' alignment positions.", sep="")

# __________________________________________________________________________________________________
if args.verbose:
	print("  Defining options...")
if args.outFile is None:
	outFile = re.sub("\\.[^\\.]+$", "_positionInfo.tsv", args.fastaFile)
else:
	outFile = args.outFile

if args.fields is None:
	fields = ['position', 'shannon', 'richness', 'unique', 'repetitions', 'cover', 'smooth']
else:
	fields = args.fields
	fields = fields.split('+')

if args.base is math.e:
	base = math.e
else:
	base = args.base

if args.rangeSmooth != 50:
	rangeSmooth = args.rangeSmooth
else:
	rangeSmooth = 50

def shannon(values):
	uniques = set(values)
	length = len(values)
	entropies = []
	for value in uniques:
		u = values.count(value)
		pi = u/float(length)
		lpi = math.log(pi, base)
		si = pi * lpi
		entropies.append(si)
	#S = -(sum(entropies))
	s = sum(entropies)
	if s == 0.0:
		S = s
	else:
		S = -(s)
	return S

# __________________________________________________________________________________________________
if args.verbose:
	print("  Calculating...")
out = {}
for field in fields:
	out[field] = list()

if args.clean is not None:
	out['shannon_clean'] = list()
	fields.append('shannon_clean')

if 'smooth' in fields:
	import statistics as st
	shan = list()

for key, value in nucleotides.items():
	if args.verbose:
		print("\r    ", key, "/", length, sep="", end="")
	out['position'].append(str(key))
	if 'shannon' in fields:
		out['shannon'].append(str(shannon(value)))
	if 'richness' in fields:
		out['richness'].append(str(len(set(value))))
	if 'unique' in fields:
		out['unique'].append('|'.join(set(value)))
	if 'repetitions' in fields:
		rep = list()
		for u in set(value):
			r = str(value.count(u))
			rep.append(r)
		out['repetitions'].append('|'.join(rep))
	if 'cover' in fields:
		gaps = value.count('-')
		out['cover'].append(str(float(1 - (gaps/seqs))))
	if 'smooth' in fields:
		shan.append(shannon(value))
	if 'shannon_clean' in fields:
		values_clean = list()
		for v in value:
			if v != '-':
				values_clean.append(v)
		out['shannon_clean'].append(str(shannon(values_clean)))

# __________________________________________________________________________________________________
if 'smooth' in fields:
	if args.verbose:
		print("\n  Smoothing...")
	if (rangeSmooth * 2) > (len(shan) * 0.1):		
		print("    Warning!: You have chosen a running mean window above the ", int(((rangeSmooth*2) / len(shan))*100), "% of the length of the alignment:\n              ", len(shan), " alignment positions and a ", (rangeSmooth*2), " bp window.", sep="")
	for i in range(1, len(shan)+1):
		if i-rangeSmooth < 1:
			b = 0
		else:
			b = i-rangeSmooth
		if i+rangeSmooth > len(shan):
			e = len(shan)
		else:
			e = i+rangeSmooth
		out['smooth'].append(str(st.mean(shan[b : e])))

# __________________________________________________________________________________________________
if args.verbose:
	print("  Writing...")
with open(outFile, 'w') as outfile:
	for i in out:
		outfile.write(i + "\t")
		#print("\t".join(i), file=outfile)
		numSubItems = len(out[i])
	outfile.write("\n")
	for level in range(numSubItems):
		for i in out:
			outfile.write(str(out[i][level]) + "\t")
			#print("\t".join(str(out[i][level])) + "\t", file=outfile)
		outfile.write("\n")

# __________________________________________________________________________________________________
if args.verbose:
	print("Done")

