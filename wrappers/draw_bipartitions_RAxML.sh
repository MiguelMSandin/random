#!/bin/bash

FASTA=$1
BS_ALL=${FASTA/.fasta/_cat_BS101.RAxML_bootstrap.tre}
# BS_ALL=${FASTA/.fasta/_cat_constraintOutgroup_BS101.RAxML_bootstrap.tre}
OUT=${FASTA/.fasta/_cat_BS101}
# OUT=${FASTA/.fasta/_cat_constraintOutgroup_BS101}

echo "Concatenating"
for BS in $(ls | grep "RAxML_bootstrap\.${FASTA/.fasta/}")
do
    cat $BS >> $BS_ALL
done

echo "Searching best tree"
for INFO in $(ls | grep "RAxML_info\.${FASTA/.fasta/}")
do
    like=$(grep "Final ML Optimization Likelihood: " $INFO | sed -e 's/.* //g')
    echo "$INFO: $like" >> $OUT.likelihoods
    echo "$INFO: $like"
done

MIN=$(cut -f2 -d ":" $OUT.likelihoods | sort -n | tail -1)
BEST_TREE=$(grep "$MIN" $OUT.likelihoods | sed -e 's/: .*//g')
BEST_TREE=${BEST_TREE/RAxML_info/RAxML_bestTree}

echo "Drawing bipartitions"
# module load bioinfo-tools raxml/8.2.12-gcc
# raxmlHPC-PTHREADS-SSE3 -m GTRCAT -p $(date +%s) -f b -t $BEST_TREE -z $BS_ALL -n $OUT
raxmlHPC-PTHREADS-AVX -m GTRCAT -p $(date +%s) -f b -t $BEST_TREE -z $BS_ALL -n $OUT

echo ""
echo "Best tree found in: $BEST_TREE"

echo "" >> $OUT.likelihoods
echo "Best tree found in: $BEST_TREE" >> $OUT.likelihoods


echo "____________________________________________________________"
for i in $(ls | grep ^RAxML.*$OUT)
do
    TYPE=$(echo $i | sed -e "s/\..*//g")
    FILE=$(echo $i | sed -e "s/^[^\.]*\.//g")
    NEW=$FILE"."$TYPE
    if [[ $(echo $NEW | grep "_info") ]]
    then
        mv $i $NEW".txt"
        echo $i renamed to: $NEW".txt"
    else
        mv $i $NEW".tre"
        echo $i renamed to: $NEW".tre"
    fi
done

