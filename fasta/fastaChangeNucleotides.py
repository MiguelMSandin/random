#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Replaces nucleotides characters, such as Us to Ts or ambiguities to Ns.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input fasta file.")

parser.add_argument("-o", "--output", dest="file_out", required=False,
                    help="Output fasta file. By default will add '_changed' before the extension.")

parser.add_argument("-c", "--change", dest="change", required=False, nargs=2, default=None,
                    help="Replaces the first given character by the second given character (e.g., '-c U T' will replace all Us by Ts).")

parser.add_argument("-a", "--ambiguities", dest="ambiguities", required=False, default=None,
                    help="A shortcut to replace all ambiguities by N.")

parser.add_argument("-u", "--upper", dest="upper", required=False, action="store_true",
                    help="Outputs the nucleotides in upper cases.")

parser.add_argument("-l", "--lower", dest="lower", required=False, action="store_true",
                    help="Outputs the nucleotides in lower cases.")

args = parser.parse_args()

if args.file_out is None:
	import re
	outFile = re.sub("\\.[^\\.]+$", "_changed.", args.file_in) + re.sub(".*\\.", "", args.file_in)
else:
	outFile = args.file_out

if args.upper and args.lower:
	print("Error! Both -u/--upper and -l/--lower arguments cannot be provided. Please choose one.")
	import sys
	sys.exit(1)

with open(outFile, "w") as outfile:
	for line in open(args.file_in):
		if ">" in line:
			print(line, end="", file=outfile)
		else:
			lineout = list()
			for l in list(line):
				o = l
				if args.change is not None:
					if l == args.change[0]:
						o = args.change[1]
				if args.ambiguities is not None:
					if l == "A" or l == "C" or l == "G" or l == "T" or l == "U":
						o = l
					else:
						o = "N"
				lineout.append(o)
			lineout = "".join(lineout)
			if args.upper:
				lineout = lineout.upper()
			if args.lower:
				lineout = lineout.lower()
			print(lineout, end="", file=outfile)

