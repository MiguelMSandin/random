#!/usr/bin/env python3

import argparse
from Bio import AlignIO
import re
import sys
import pandas as pd

parser = argparse.ArgumentParser(description="Builds the consensus sequence of an alignment.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-f", "--file", dest="inFile", required=True,
                    help="An aligned fasta file.")

parser.add_argument("-o", "--output", dest="outFile", required=False, default=None,
                    help="The output name of the file to write the consensus sequence (in fasta format). If the file exists, it will append the consensus sequence at the end of the file. By default, the consensus sequence will be printed in the console.")

parser.add_argument("-t", "--threshold", dest="threshold", required=False, action='store', type=int, default=70,
                    help="The threshold to build a consensus given in percentage (0-100). Default=70.")

parser.add_argument("-b", "--base", dest="baseThreshold", required=False, action='store', type=int, default=30,
                    help="The threshold for a base to be considered. Default=30")

parser.add_argument("-g", "--gaps", dest="gaps", required=False, action='store', type=int, default=80,
                    help="If the proportion of gaps is more than 'g', the consensus is a gap ('-'). Default=80")

parser.add_argument("-a", "--ambiguities", dest="ambiguities", required=False, default=None, action="store_true",
                    help="If selected, ambiguites will be set to 'N' and not the IUPAC bases.")

parser.add_argument("-m", "--most", dest="most", required=False, default=None, action="store_true",
                    help="If selected, the most abundant sequence at each position will be used instead of the IUPAC code.")

parser.add_argument("-r", "--removeGaps", dest="removeGaps", required=False, default=None, action="store_true",
                    help="If selected, gaps in the consensus sequence will be remove.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, default=None, action="store_true",
                    help="If selected, will print information in the console.")

args = parser.parse_args()

# __________________________________________________________________________________________________
if args.verbose:
	print("  Setting variables")
	print("    Alignment:  ", args.inFile)
	print("    Threshold:   0.", args.threshold, sep="")
	print("    Gaps:        0.", args.gaps, sep="")
	print("    Base:        0.", args.baseThreshold, sep="")
	if args.outFile is not None:
		print("    Output file: ", args.outFile, sep="")
	if args.ambiguities is not None:
		print("    Ambiguities will be set to 'N'")
	if args.most is not None:
		print("    The most abundant base will be taken for consensus")
	if args.removeGaps is not None:
		print("    Gaps will be removed in output sequence")

iupac = {"A":"A", "C":"C", "G":"G", "T":"T", "AG":"R", "GA":"R", "CT":"Y", "TC":"Y", "GC":"S", "CG":"S", "AT":"W", "TA":"W", "GT":"K", "TG":"K", "AC":"M", "CA":"M", "CGT":"B", "CTG":"B", "TCG":"B", "GCT":"B", "TGC":"B", "GTC":"B", "AGT":"D", "ATG":"D", "GAT":"D", "TAG":"D", "GTA":"D", "TGA":"D", "ACT":"H", "ATC":"H", "CAT":"H", "TAC":"H", "CTA":"H", "TCA":"H", "ACG":"V", "AGC":"V", "CAG":"V", "GAC":"V", "CGA":"V", "GCA":"V", "ACTG":"N", "ACGT":"N", "ATCG":"N", "AGCT":"N", "ATGC":"N", "AGTC":"N", "CATG":"N", "CAGT":"N", "TACG":"N", "GACT":"N", "TAGC":"N", "GATC":"N", "CTAG":"N", "CGAT":"N", "TCAG":"N", "GCAT":"N", "TGAC":"N", "GTAC":"N", "CTGA":"N", "CGTA":"N", "TCGA":"N", "GCTA":"N", "TGCA":"N", "GTCA":"N"}

# __________________________________________________________________________________________________
if args.verbose:
	print("  Reading alignment...", end="")

fasta =  AlignIO.read(args.inFile, "fasta")

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
	print("\r  Reading alignment and asigning positions...", end="")

nucleotides = {}
for position in range(length):
	p = position +1
	nucleotides[p] = list(fasta[:, position].upper())

if args.verbose:
	print("\r    Sequences in alignment:", seqs, "                    ")
	print("    Positions in alignment:", length)

# __________________________________________________________________________________________________
if args.verbose:
	print("  Calculating...")

consensus = {}
consensus['consensus'] = list()
if args.ambiguities is not None:
	consensus['ambiguities'] = list()
if args.gaps is not None:
	consensus['gaps'] = list()
if args.ambiguities is not None and args.gaps is not None:
	consensus['consensusGaps'] = list()
if args.most is not None:
	consensus['most'] = list()

for position, value in nucleotides.items():
	if args.verbose:
		print("\r    ", position, "/", length, sep="", end="")
	
	# Create a dataframe with the bases and its frequencies in each position
	bases = {}
	for u in set(value):
		r = int(value.count(u))
		bases[u] = r
	bases = pd.DataFrame({'bases':list(bases.keys()),'occurrence':list(bases.values())})
	bases['percentage'] = (bases['occurrence']) / sum(bases['occurrence'])
	
	# First dealing with the gaps
	if '-' in list(bases['bases']):
		if float(bases[bases['bases']=='-']['percentage']) >= args.gaps/100:
			base = list('-')
		else: 
			bases = bases[bases['bases']!='-']
			bases['percentage'] = (bases['occurrence']) / sum(bases['occurrence'])
	
	# if there is only one base, take its value in the consensus
	if len(bases['bases']) == 1:
		base = bases['bases']
	
	# dealing now with the consensus
	if args.most is not None:
		baseMost = bases['bases'][bases['percentage']==max(bases['percentage'])]
		if len(baseMost) != 1:
			baseMost=iupac[''.join(baseMost)]
		consensus['most'] = (consensus['most'] + list(baseMost))
	if len(bases['bases']) > 1:
		bases = bases[bases['percentage'] > args.baseThreshold/100]
		bases['percentage'] = (bases['occurrence']) / sum(bases['occurrence'])
		if len(list(bases['bases'][bases['percentage'] >= args.threshold/100])) == 0:
			b = list(bases['bases'])
			b = ''.join(b)
			for k, v in iupac.items():
				if b == k:
					base = v
		elif len(list(bases['bases'][bases['percentage'] >= args.threshold/100])) == 1:
			base = list(bases['bases'][bases['percentage'] >= args.threshold/100])
		elif len(list(bases['bases'][bases['percentage'] >= args.threshold/100])) > 1:
			b = list(bases['bases'][bases['percentage'] >= args.threshold/100])
			b = ''.join(b)
			for k, v in iupac.items():
				if b == k:
					base = v
	consensus['consensus'] = (consensus['consensus'] + list(base))

if args.verbose:
	print("\n")

consensus['consensus'] = ''.join(consensus['consensus'])
ext=""
if args.ambiguities is not None:
	consensus['ambiguities'] = re.sub("[^ACTG-]", "N", consensus['consensus'])
if args.gaps is not None:
	consensus['gaps'] = re.sub("-", "", consensus['consensus'])
if args.ambiguities is not None and args.gaps is not None:
	consensus['consensusGaps'] = re.sub("-", "", consensus['consensus'])
	consensus['consensusGaps'] = re.sub("[^ACTG]", "N", consensus['consensusGaps'])
if args.most is not None:
	consensus['most'] = ''.join(consensus['most'])
	ext = "_mostAbundant"

# __________________________________________________________________________________________________

if args.ambiguities is None and args.removeGaps is None and args.most is None:
	sequence = str(consensus['consensus'])
elif args.most is not None:
	if args.removeGaps is not None:
		sequence = str(re.sub("-", "", consensus['most']))
	else:
		sequence = str(consensus['most'])
elif args.ambiguities is not None and args.removeGaps is not None:
	sequence = str(consensus['consensusGaps'])
elif args.ambiguities is not None:
	sequence = str(consensus['ambiguities'])
elif args.removeGaps is not None:
	sequence = str(consensus['gaps'])

if args.outFile is not None:
	f = open(args.outFile, "a")
	tmp = re.sub("\\.[^\\.]+$", "", args.inFile)
	f.write(str(">" + str(tmp) + "_consensus_t" + str(args.threshold) + "_b" + str(args.baseThreshold) + "_g" + str(args.gaps) + ext + "\n"))
	f.write(str(str(sequence) + "\n"))
	f.close()
else:
	print("  Consensus sequence:")
	print("")
	print(str(sequence))
	print("")

if args.verbose:
	print("  Consensus positions:   ", len(sequence)-sequence.count("-"))
	amb = len(sequence)-len(re.sub("[^ACTG-]", "", sequence))
	if amb > 0:
		print("  Of which are ambiguous:", amb)
	if args.removeGaps is None:
		print("  Of which are gaps:     ", sequence.count("-"))
	print("")
	print("Done")

