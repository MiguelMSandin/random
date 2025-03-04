#!/usr/bin/env python3

import argparse
import re

version='0.4.0-beta'

parser = argparse.ArgumentParser(description="%s v%s: Select sequence in a fasta file either from a list or that match a pattern and extract them or remove them." % ('%(prog)s', version), add_help=False)

requirArgs = parser.add_argument_group('Required arguments')
eitherArgs = parser.add_argument_group('Required at least one of')
outputArgs = parser.add_argument_group('Optional arguments related to the output')
optionArgs = parser.add_argument_group('Other optional arguments')

requirArgs.add_argument("-f", "--file", dest="file_in", required=True,
						help="Input fasta file. Avoid spaces in the sequence names.")

eitherArgs.add_argument("-l", "--list", dest="listSeq", required=False, default=None,
						help="List of sequences to be selected. This must be a different file with each sequence name in a different line.")

eitherArgs.add_argument("-p", "--pattern", dest="pattern", required=False, default=None, nargs="+",
						help="Pattern(s) to be matched for selection of the sequences. When using this option in combination with '-k/--keep', it might be faster to use 'grep -A 1 PATTERN FILE_IN > FILE_OUT' if each sequence from the input file is in a single line.")

outputArgs.add_argument("-o", "--output", dest="file_out", required=False, default=None,
						help="Output file. By default will add '_selected' to the input file name. If the file already exists, sequences will be appended at the end of the file.")

outputArgs.add_argument("-k", "--keep", dest="keep", required=False, action="store_true",
						help="If selected, will export selected sequences. If neither this argument nor '-r/--remove' are selected, this argument will be set by default.")

outputArgs.add_argument("-r", "--remove", dest="remove", required=False, action="store_true",
						help="If selected, will delete selected sequences. If both this argument and '-k/--keep' are selected, this argument will be ignored.")

optionArgs.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
						help="If selected, will not print information to the console.")

optionArgs.add_argument("-h", "--help", action="help",
						help="Show this help message and exit.")

optionArgs.add_argument("-V", "--version", action='version',
						version='oligoN-design %s v%s' % ('%(prog)s', version),
						help='Show the version number and exit.')

args = parser.parse_args()

# Define functions _________________________________________________________________________________
def readFasta(fastafile):
	out = {}
	for line in open(fastafile):
		if line.startswith(">"):
			name = line.replace(">", "")
			name = name.replace("\n", "")
			out[name] = str()
		else:
			sequence = line.replace("\n", "")
			sequence = sequence.replace("-", "")
			sequence = sequence.upper()
			out[name] = (out[name] + sequence)
	return out

# Troubleshoot input variables _____________________________________________________________________
if args.listSeq is None and args.pattern is None:
	import sys
	print("Error: You should specify either a list (-l/--list) or a pattern (-p/--pattern) so I can select sequences...")
	sys.exit(1)

keep = False
remove = False
if args.keep is False and args.remove is False:
	keep = True
elif args.keep:
	keep = True
elif args.remove:
	remove = True

if args.keep and args.remove and args.verbose:
	print("  Warning! '-r/--remove' has been disabled")

# Output file
if args.file_out is None:
	outFile = re.sub("\\.[^\\.]+$", "_selected", args.file_in) + re.sub(".*\\.", ".", args.file_in)
else:
	outFile = args.file_out

# Reading input files ______________________________________________________________________________
if args.listSeq is not None:
	listSeq = [line.strip() for line in open(args.listSeq)]

if args.verbose:
	print("  Reading input file:    ", args.file_in)
infile = readFasta(args.file_in)

# Selecting sequences ______________________________________________________________________________
seq_in = 0
seq_out = 0
with open(outFile, "a") as outfile:
	seqsid = set()
	for key, val in infile.items():
		seq_in += 1
		seqsid.add(key)
		if keep:
			if args.listSeq is not None:
				if key in listSeq:
					seq_out += 1
					print(">" + str(key) + "\n" + str(val), file=outfile, flush=True)
			if args.pattern is not None:
				for pattern in args.pattern:
					if pattern in key:
						seq_out += 1
						print(">" + str(key) + "\n" + str(val), file=outfile, flush=True)
		elif remove:
			if args.listSeq is not None:
				if key not in listSeq:
					seq_out += 1
					print(">" + str(key) + "\n" + str(val), file=outfile, flush=True)
			if args.pattern is not None:
				out = True
				for pattern in args.pattern:
					if pattern in key:
						out = False
				if out:
					seq_out += 1
					print(">" + str(key) + "\n" + str(val), file=outfile, flush=True)

if args.verbose:
	print("  Output file written to:", outFile)
	if args.listSeq is not None:
		missing = set()
		for l in listSeq:
			if l not in seqsid:
				missing.add(l)
		if len(missing) > 0:
			print("  The following sequences from the list are not found in the input fasta file:\n", "\n ".join(missing))
	print("    Sequences in: ", seq_in)
	print("    Sequences out:", seq_out)
	print("Done")
