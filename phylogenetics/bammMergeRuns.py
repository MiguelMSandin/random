#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Combines several BAMM runs into single files. This script assumes the files to be combined are correct, since it will not check for errors (e.g.; different columns or different sampling frequency).")

requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--files", dest="files", required=True, nargs='+',
                    help="Input files to be concatenated in the given order.")

requiredArgs.add_argument("-o", "--output", dest="fileOut", required=True,
                    help="Output file.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

i = 0
generation = 0
generationh = -1
step = 0
with open(args.fileOut, "w") as outfile:
	for filei in args.files:
		i += 1
		j = 0
		pct = -1
		lines = sum(1 for line in open(filei))
		for line in open(filei):
			j += 1
			if args.verbose:
				pcti = round(j/lines*100)
				if pcti > pct:
					pct = pcti
					print("\r  Concatenating file ", filei, "\t(", i, "/", len(args.files), ")\t", pct, "%", sep="", end="", flush=True)
			if i == 1 and j == 1:
				header =  line.strip().split(",")
				print(','.join(header), file=outfile, flush=True)
			elif i == 1 and j !=1:
				lineSplit =  line.strip().split(",")
				print(','.join(lineSplit), file=outfile, flush=True)
				generationi = int(lineSplit[0])
				if generationi > generation:
					if generationi - generation > step:
						step = generationi - generation
					generation = generationi
			elif i > 1 and j != 1:
				lineSplit = line.strip().split(",")
				generationi = int(lineSplit[0])
				if generationi > generation:
					if generationi - generation > step:
						step = generationi - generation
					generation = generationi
				if generationi < generation and generationh == generationi:
					generation = generation
				if generationi < generation and generationh != generationi:
					generation = generation + step
				lineout = list()
				lineout.append(str(generation))
				lineSplit.pop(0)
				lineout.extend(lineSplit)
				lineout = ','.join(lineout)
				print(lineout, file=outfile, flush=True)
				generationh = generationi
		if args.verbose:
			print("")

if args.verbose:
	print("Done")
