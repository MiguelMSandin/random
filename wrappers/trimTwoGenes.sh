#!/bin/bash

usage()
{
    echo ""
    echo "A wrapper to trim independently two concatenated genes within a fasta file."
    echo ""
    echo "Requires the scripts 'fastaSplit.py', 'fastaConcat.py' and the package 'trimAl'" 
    echo "(from Capella-Gutierrez Silla-Martinez, Gabaldon. Bioinformatics 2009, 25:1972-1973.)"    
    echo ""
    echo "Usage: trimTwoGenes.sh -f fastaFile -p position -g gapThreshold"
    echo ""
    echo "  -h    Print this information"
    echo ""
    echo "  -f    A fasta file where two genes have been concatenated"
    echo "  -p    The position where the first gene ends. As in 'fastaSplit.py -p'"
    echo "  -g    1 - (fraction of sequences with a gap allowed). As in 'trimal -gt'"
    echo ""
}

while getopts "hf:p:g:" opt; do
	case ${opt} in
		h )
			usage
			exit
			;;
		f )
			FASTA=$OPTARG
			;;
		p )
			POSITION=$OPTARG
			;;
		g )
			GAPTHRESHOLD=$OPTARG
			;;
	esac
done

echo "Splittingfile"
fastaSplit.py -f $FASTA -p "+$POSITION+" -o tmp -r

echo "Trimming"
for FILE in $(ls tmp)
do
	trimal -in tmp/$FILE -out tmp/${FILE/.fasta/_trim.fasta} -gt $GAPTHRESHOLD
	rm tmp/$FILE
done

echo "Concatenating"
fastaConcat.py -f tmp/${FASTA/.fasta/_1_trim.fasta} tmp/${FASTA/.fasta/_2_trim.fasta} -o ${FASTA/.fasta/_trim.fasta} -a

rm -fr tmp
