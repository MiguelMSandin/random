#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import re
import statistics
import numpy as np

parser = argparse.ArgumentParser(description="Returns overall statistics and numbers from a fasta file.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="fileIn", nargs='+', required=True,
                    help="An aligned fasta file.")

parser.add_argument("-s", "--short", dest="short", required=False, action="store_true",
					help="If selected, will only print number of sequences, total number of bases (and if aligned) aligned positions and proportion of gaps.")

parser.add_argument("-d", "--detailed", dest="detailed", required=False, action="store_true",
					help="If selected, besides printing the overall statistics, will print similar statistics for every sequence. This option is incopatible with the '-s/--short' option.")

args = parser.parse_args()

if(args.short):
	print("File\tSequences\tBases\t(Positions\tGaps)")

for file in args.fileIn:
	seqs = 0
	lengths = []
	lengthsRaw = []
	As = 0
	Cs = 0
	Gs = 0
	Ts = 0
	Ns = 0
	gaps = 0
	ambiguities = list()
	for i in SeqIO.parse(open(file), "fasta"):
		seqs += 1
		seq = i.seq
		lengthsRaw.append(len(seq))
		lengths.append(len(re.sub("-", "", str(seq))))
		a = seq.upper().count("A")
		As += a
		c = seq.upper().count("C")
		Cs += c
		g = seq.upper().count("G")
		Gs += g
		t = seq.upper().count("T")
		Ts += t
		gap = seq.upper().count("-")
		gaps += gap
		n = len(seq) - (a + c + g + t + gap) # Count number of ambiguities
		Ns += n
		if n > 0:
			tmp = re.sub("A|C|G|T|-", "", str(seq.upper()))
			ambiguities.append(list(tmp))
	
	ambiguities = [item for l in ambiguities for item in l]
	ambiguities = set(ambiguities)
	
	totalbp = As + Cs + Gs + Ts + Ns
	
	if(args.short):
		if(len(set(lengthsRaw)) > 1):
			print(str(file), ":\t", seqs, "\t", totalbp, sep="")
		else:
			totalPositions = int(*set(lengthsRaw)) * seqs
			gaps = totalPositions - totalbp
			print(str(file), ":\t", seqs, "\t", totalbp, "\t", *set(lengthsRaw), "\t", round(gaps / totalPositions * 100, 2), "%", sep="")
	else:
		print("File:\t", file)
		if(len(set(lengthsRaw)) > 1):
			print("  Not aligned")
		else:
			print("  Aligned positions:\t", *set(lengthsRaw))
		print("  Number of sequences:\t", seqs)
		print("  Shortest sequence:   \t", min(lengths), "bp")
		print("    5th percentile:    \t", int(np.percentile(lengths, 5)), "bp")
		print("    25th percentile:   \t", int(np.percentile(lengths, 25)), "bp")
		print("    Median:            \t", int(np.percentile(lengths, 50)), "bp")
		print("  Average length:      \t", round(statistics.mean(lengths), 2), "bp")
		print("    Standard deviation:\t", round(statistics.stdev(lengths), 2), "bp")
		print("    75th percentile:   \t", int(np.percentile(lengths, 75)), "bp")
		print("    95th percentile:   \t", int(np.percentile(lengths, 95)), "bp")
		print("  Longest sequence:    \t", max(lengths), "bp")
		print("  Base composition:")
		print("    A:\t", As, "bp\t", round(As/totalbp*100, 2), "%")
		print("    C:\t", Cs, "bp\t", round(Cs/totalbp*100, 2), "%")
		print("    G:\t", Gs, "bp\t", round(Gs/totalbp*100, 2), "%")
		print("    T:\t", Ts, "bp\t", round(Ts/totalbp*100, 2), "%")
		if Ns > 0:
			print("    Ambiguities:\t", Ns, "bp\t", round(Ns/totalbp*100, 2), "%. Bases:", ", ".join(ambiguities))
		else:
			print("    No ambiguities found")
		if(len(set(lengthsRaw)) > 1):
			print("  Total bases:\t", totalbp)
		else:
			totalPositions = int(*set(lengthsRaw)) * seqs
			gaps = totalPositions - totalbp
			print("  Gaps in alignment:\t", gaps, "\t", round(gaps / totalPositions * 100, 2), "%")
			print("  Total bases:      \t", totalbp, "\t", round( totalbp / totalPositions * 100, 2), "%")
		print("")
		if args.detailed:
			if(len(set(lengthsRaw)) > 1):
				print("sequence\tlength\tA\tC\tG\tT\tambiguities")
			else:
				print("sequence\tlength\tA\tC\tG\tT\tambiguities\tgaps\tproportionGaps")
			for i in SeqIO.parse(open(file), "fasta"):
				name = i.id
				seq = i.seq
				l = len(re.sub("-", "", str(seq)))
				a = seq.upper().count("A")
				c = seq.upper().count("C")
				g = seq.upper().count("G")
				t = seq.upper().count("T")
				gap = seq.upper().count("-")
				n = len(seq) - (a + c + g + t + gap) # Count number of ambiguities
				gapp = round(gap / len(seq) * 100, 2)
				if(len(set(lengthsRaw)) > 1):
					print(str(name) + '\t' + str(l) + '\t' + str(a) + '\t' + str(c) + '\t' + str(g) + '\t' + str(t) + '\t' + str(n))
				else:
					print(str(name) + '\t' + str(l) + '\t' + str(a) + '\t' + str(c) + '\t' + str(g) + '\t' + str(t) + '\t' + str(n) + '\t' + str(gap) + '\t' + str(gapp))

