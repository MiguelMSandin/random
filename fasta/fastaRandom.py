#!/usr/bin/env python3

import argparse
import random
import sys

parser = argparse.ArgumentParser(description="Generates a random fasta file of s sequences of length l.")

# Add the arguments to the parser
parser.add_argument("-s", "--sequences", dest="sequences", required=False, default=1, type=int,
                    help="The number of sequences to output. Default = 1")

parser.add_argument("-l", "--length", dest="length", required=False, default=1000, type=int,
					help="The length of the sequences in base pairs. Default = 1000")

parser.add_argument("-o", "--output", dest="file_out", required=False, default=None,
                    help="Output file, if exists, output will be concatenated at the end. If not provided, it will print the fasta file to the console.")

parser.add_argument("-n", "--normal", dest="normal", required=False, default=None,
					help="If selected, will adjust the sequence length to a normal distribution of length l ('-l/--length') and standard deviation parsed in this command.")

parser.add_argument("-m", "--mutations", dest="mutations", required=False, default=None,
					help="If selected, will add random mutations to the first randomly generated sequence with the given probability (i.e.; -m 0.1)")

parser.add_argument("-p", "--protein", dest="protein", required=False, action="store_false",
					help="If selected, will generate a protein sequence instead.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_true",
					help="If selected, will print information to the console.")

args = parser.parse_args()

if args.protein:
	choices = ["A", "C", "G", "T"]
	sequence = "DNA"
else:
	choices = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "U", "O"]
	sequence = "protein"

if args.normal is not None and args.mutations is not None:
	print("\nAt the moment, a normal distribution and adding mutations are not compatible.\nPlease remove one of the two options.\n")
	sys.exit(1)

if args.normal is not None:
	import numpy as np
	rng = np.random.default_rng()

if args.verbose and args.file_out is not None:
	if args.normal is None and args.mutations is None:
		print("  Exporting", args.sequences, "random", sequence, "sequences of length", args.length, "to", args.file_out)
	if args.normal is not None:
		print("  Exporting ", args.sequences, " random ", sequence, " sequences of length N(", args.length, ", ", args.normal, ") to ", args.file_out, sep="")
	if args.mutations is not None:
		print("  Exporting", args.sequences, sequence, "sequences of length", args.length, "and a", args.mutations, "probability of mutations to", args.file_out)

if args.mutations is not None:
	if args.sequences == 1:
		print("Warning! Only one sequence has been selected, there will be no more sequences to add mutations to...")
	mutationRate = float(args.mutations)
	template = random.choices(choices, k=args.length)
	out = "".join(template)
	if args.file_out is not None:
		f = open(args.file_out, "a")
		f.write(str(">r1\n"))
		f.write(str(str(out) + "\n"))
		f.close()
	else:
		print(">r1", sep="")
		print(str(out))
	for i in range(args.sequences-1):
		if args.verbose and args.file_out is not None:
			print("\r    ", round((i+2)/args.sequences*100), "%", sep="", end="")
		out = list()
		for j in template:
			mutation = random.choices([1, 0], weights=[mutationRate, 1-mutationRate])
			if mutation == [1]:
				m = random.choices(choices)
				out = out + m
			else:
				out.append(j)
		out = "".join(out)
		name = str("r") + str(i+2)
		if args.file_out is not None:
			f = open(args.file_out, "a")
			f.write(str(">" + name + "\n"))
			f.write(str(str(out) + "\n"))
			f.close()
		else:
			print(">", name, sep="")
			print(str(out))
else:
	for i in range(args.sequences):
		if args.verbose and args.file_out is not None:
			print("\r    ", round((i+1)/args.sequences*100), "%", sep="", end="")
		if args.normal is None:
			l = args.length
		else:
			l = round(rng.normal(args.length, float(args.normal)))
		out = random.choices(choices, k=l)
		out = "".join(out)
		name = str("r") + str(i+1)
		if args.file_out is not None:
			f = open(args.file_out, "a")
			f.write(str(">" + name + "\n"))
			f.write(str(str(out) + "\n"))
			f.close()
		else:
			print(">", name, sep="")
			print(str(out))

if args.verbose:
	print("\n", sep="", end="")
