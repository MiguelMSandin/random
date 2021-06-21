#!/usr/bin/env python3

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="From a targeted fasta file, select all possible primers of length 'l', that matches 'm' sequences and has a maximum specificity 's' to the reference database.")

# Add the arguments to the parser
parser.add_argument("-t", "--target", dest="target", required=True,
					help="A fasta file with the sequences you want to find a primer to match.")

parser.add_argument("-r", "--reference", dest="reference", required=True,
					help="A reference file to look against. The targeted group shouldn't be included in the reference. If you are not sure whether it is included or not use the option '-s/--specificity'.")

parser.add_argument("-o", "--output", dest="output", required=True,
					help="The name of the output fasta and log files. Please note the extensions '.fasta' and '.tsv' will be added to the specify name respectively.")

parser.add_argument("-l", "--length", dest="length", required=False, action='store', default='18+22',
					help="The desire length of the primers to be searched. Several lengths can be selected by specifying a range with a '+' sign in between. By default it will look for primers of length 18, 19, 20, 21 and 22 base pairs ('18+22').")

parser.add_argument("-s", "--specificity", dest="specificity", required=False, action='store', type=float, default=0.01,
					help="The maximum percentage of sequences that can be matched in the reference file (0 >= s >= 1). Default = 0.01")

parser.add_argument("-S", "--specificityAbs", dest="specificityAbs", required=False, action='store', type=int, default=None,
					help="Same as '-s/--specificity' but absolute values.")

parser.add_argument("-m", "--minimum", dest="minimum", required=False, action='store', type=float, default=0.8,
					help="The minimum percentage of sequences that the primer has to match in the target file (0 >= m >= 1). Default = 0.8")

parser.add_argument("-M", "--minimumAbs", dest="minimumAbs", required=False, action='store', type=int, default=None,
					help="Same as '-m/--minimum' but absolute values.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, default=None, action="store_true",
					help="If selected, will print information to the console.")

args = parser.parse_args()

if args.verbose is not None:
	verbose = True
	print("  Setting variables...")
else:
	verbose = False

# Setting the range of lengths _____________________________________________________________________
if '+' in args.length:
	mi = args.length.split('+')[0]
	ma = args.length.split('+')[1]
	lengths = range(int(mi), int(ma)+1)
else:
	lengths = int(args.lengths)
if verbose:
	print("    Lengths:    ", *lengths)

# Setting the specificities ________________________________________________________________________
if args.specificityAbs is not None:
	specificity = args.specificityAbs
elif args.specificity != 0.01:
	specificity = args.specificity
else:
	specificity = 0.01
	if args.specificityAbs is not None:
		print("    Warning!! '-s/--specificity' has been parsed by '-S/--specificityAbs'. Taking absolute values.")
if verbose:
	print("    Specificity:", specificity)

if args.minimumAbs is not None:
	minimum = args.minimumAbs
elif args.minimum != 0.8:
	minimum = args.minimum
else:
	minimum = 0.8
	if args.minimumAbs is not None:
		print("    Warning!! '-m/--minimum' has been parsed by '-M/--minimumAbs'. Taking absolute values.")
if verbose:
	print("    Minimum:    ", minimum)

# Reading target file ______________________________________________________________________________
if verbose:
	print("  Reading target file...")

target = {}
w1 = 0
for line in SeqIO.parse(open(args.target), "fasta"):
	target[line.id] = str(line.seq.upper())
	if "-" in line.seq:
		w1 += 1

if verbose:
	print("    Sequences in target file:   ", len(target))
	if w1 > 0:
		print("    Warning!!", args.target, "contains gaps")

# Reading reference file ___________________________________________________________________________
if verbose:
	print("  Reading reference file...")

ref = {}
w2 = 0
for line in SeqIO.parse(open(args.reference), "fasta"):
	ref[line.id] = str(line.seq.upper())
	if "-" in line.seq:
		w2 += 1

if verbose:
	print("    Sequences in reference file:", len(ref))
	if w2 > 0:
		print("    Warning!!", args.reference, "contains gaps")

# Start the search _________________________________________________________________________________
primers = set()
logFile = str(args.output + ".tsv")
fasFile = str(args.output + ".fasta")
with open(logFile, "w") as logfile, open(fasFile, "w") as fasfile:
	print("identifier\tsequence\tmatch_ref\tmatch_ref_abs\tmatch_target\tmatch_target_abs\t\t", file=logfile)
	for length in lengths:  # Loop through the different lengths
		if verbose:
			print("  Searching primers of", str(length), "base pairs...")
			i = 0
		for tv in target.values():  # Loop through all the sequences in the target file
			if verbose:
				i += 1
				print("\r    ", str(round(i/len(target)*100, 1)), "%", sep="", end="")
				if i == len(target):
					print("\r    Completed")
			l = len(tv)
			for p in range(0, l-length+1):  # Split the target ith sequence into potential primers of length 'l' 
				pprimer = tv[int(p):(int(p+length))]
				if pprimer not in primers:  # Check if the potential primer is already in the list
					r = 0
					for toCheck in target.values():  # Check how many times the potential primer is repeated in the target file
						if pprimer in toCheck:
							r += 1
					if args.minimumAbs is not None:  # Check if the potential primer is repeated more than 'M' times
						R = r
					else:  # Check if the potential primer is repeated more than 'm' times
						R = r/len(target)
					if R >= minimum:  # Check if the potential primer matches the requirements against the target file
						c = 0
						for toCheck in ref.values():  # Check how many times the potential primer is repeated in the reference file
							if pprimer in toCheck:
								c += 1
					if args.specificityAbs is not None:  # Check if the potential primer is repeated less than 'S' times
						C = c
					else:  # Check if the potential primer is repeated less than 's' times
						C = c/len(ref)
					if C <= specificity:  # Check if the potential primer matches the requirements against the reference file
							primers.add(pprimer)
							print("primer" + str(len(primers)) + "\t" + str(pprimer) + "\t" + str(c/len(ref)) + "\t" + str(c) + "\t" + str(r/len(target)) + "\t" + str(r), file=logfile)
							print(">primer" + str(len(primers)) + "\n" + pprimer, file=fasfile)


