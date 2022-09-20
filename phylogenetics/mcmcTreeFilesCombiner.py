#!/usr/bin/env python3

import argparse
import sys

parser = argparse.ArgumentParser(description="Removes reciprocal (A-B = B-A) and identical (A-A) hits from a tab separated list. The first and second columns are taken as identifiers. The rest of the columns are irrelevant.")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--files", dest="files", required=True, nargs='+',
                    help="Input txt files to be cleaned and bind in the given order.")

requiredArgs.add_argument("-o", "--output", dest="fileOut", required=True,
                    help="Output file.")

parser.add_argument("-b", "--burnin", dest="burnin", required=False, type=int, default="10",
                    help="The percentage of samples to be removed at the beginning of each file. Default=10.")

parser.add_argument("-s", "--states", dest="states", required=False, type=int, default=None,
                    help="The number of samples to be removed at the beginning of each file. If selected, will ignore the percentage selected by the option burnin. Default=None")

parser.add_argument("-c", "--clean", dest="clean", required=False, action="store_false",
                    help="If selected, will not remove incomplete lines.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

warning = 0
state = 0
i = 0
with open(args.fileOut, "w") as outfile:
	for filei in args.files:
		i += 1
		lines = sum(1 for line in open(filei)) - 1
		if args.verbose:
			print("  Processing file ", filei, " (", i, "/", len(args.files), ")", sep="")
			if args.states is not None:
				burn = args.states
				print("    Burnin: ", burn, " states or ", round(burn / lines * 100, 2), "%", sep="")
			else:
				burn = int(lines * args.burnin / 100)
				print("    Burnin: ", args.burnin, "% or ", burn, " states", sep="")
			print("    States out:", lines - burn)
		j = 0
		for line in open(filei):
			j += 1
			if j == 1:
				if i == 1:
					print(line, end="", file=outfile)
					headers = line.strip().split()
				else:
					tmp = line.strip().split()
					if len(tmp) != len(headers):
						print("\nError: File", filei, "have different number of headers, please check all files are coming from the same analysis.\nStopping\n")
						sys.exit(1)
			elif j >= burn:
				tmp = line.strip().split()
				if args.clean:
					if len(tmp) == len(headers):
						state += 1
						tmp[0] = str(state)
						lineOut = '\t'.join(tmp)
						print(lineOut, file=outfile)
				else:
					if len(tmp) != len(headers):
						warning += 1
					state += 1
					tmp[0] = str(state)
					lineOut = '\t'.join(tmp)
					print(lineOut, file=outfile)

if warning != 0:
	print("  Warning! There are", warning, "incomplete lines in the output file")

if args.verbose:
	print("  In total", len(args.files), "files were combined into", state, "total states")
	print("Done")
