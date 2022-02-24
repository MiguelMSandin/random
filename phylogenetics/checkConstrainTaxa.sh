#!/bin/bash

usage()
{
    echo ""
    echo "Usage: checkConstrainTaxa.sh -d DIRECTORY -o OUTPUT"
    echo ""
    echo "  -d    a directory containing the selected taxa for each taxon."
    echo "  -o    an output file containing the unique summarized taxons at a higher taxonomic level. If exists, will be overwritten."
    echo ""
    echo "Please check the code according to the names of the sequences..."
    echo ""
}

while getopts "hd:o:" opt; do
        case ${opt} in
                h )
                        usage
                        exit
                        ;;
                d )
                        DIRECTORY=$OPTARG
                        ;;
                o )
                        OUT=$OPTARG
                        ;;
        esac
done

if test -f "$OUT"; then
	echo "$OUT file exists, deleting..."
	rm -f $OUT
fi

echo "Searching in $DIRECTORY"

tmp1=$(mktemp --tmpdir=$(pwd))
tmp2=$(mktemp --tmpdir=$(pwd))

i=0
c=0
LENGTH=$(ls -1 $DIRECTORY | wc -l)
for FILE in $(ls $DIRECTORY)
do
	((i=i+1))
	echo -n -e "\r  $i/$LENGTH: $FILE                    "
	
	if [ -s "$DIRECTORY/$FILE" ]; then
		
		tmpc1=$(grep -c "PR2_" "$DIRECTORY/$FILE")
		if [ $tmpc1 != 0 ]; then
			grep "PR2_" "$DIRECTORY/$FILE" | \
				cut -d'|' -f 3,4,5 | \
				sort | \
				uniq | \
				sed 's/|/\t/g' > $tmp1
# 			if [ $(wc -l < $tmp1) != 1 ]; then
# 				echo "$FILE" >> $OUT
# 				echo "	PR2 search: $tmpc1 sequences" >> $OUT
# 				cat $tmp1 >> $OUT
# 				echo "" >> $OUT
# 			fi
		fi
		
		tmpc2=$(grep -c "PacBio_" "$DIRECTORY/$FILE")
		if [ $tmpc2 != 0 ]; then
			grep "PacBio_" "$DIRECTORY/$FILE" | \
				sed 's/_X/-X/g' | \
				sed 's/.*Eukaryota//g' | \
				cut -d'_' -f 3,4,5 | \
				sort | \
				uniq | \
				sed 's/_/\t/g' > $tmp2
# 			if [ $(wc -l < $tmp2) != 1 ]; then
# 				echo "$FILE" >> $OUT
# 				echo "	PacBio search: $tmpc2 sequences" >> $OUT
# 				cat $tmp2 >> $OUT
# 				echo "" >> $OUT
# 			fi
		fi
		
		if [[ $(wc -l < $tmp1) -ge 2 || $(wc -l < $tmp2) -ge 2 ]]; then
			((c=c+1))
			echo "$FILE" >> $OUT
			if [ $(wc -l < $tmp1) -ge 2 ]; then
				echo "	PR2 search: $tmpc1 sequences" >> $OUT
				cat $tmp1 >> $OUT
			fi
			if [ $(wc -l < $tmp2) -ge 2 ]; then
				echo "	PacBio search: $tmpc2 sequences" >> $OUT
				cat $tmp2 >> $OUT
			fi
			echo "" >> $OUT
		fi
		unset tmpc1 tmpc2
		truncate -s 0 $tmp1
		truncate -s 0 $tmp2
	else
		rm -f "$DIRECTORY/$FILE"
	fi
done

echo -e "\r  Finished checking                           "

if [ c != 0 ]; then
	echo ""
	echo "  There are $c conflicting files"
	echo "  Please, check '$OUT' and correct the files in '$DIRECTORY' if neccessary"
	echo ""
else
	echo ""
	echo "  Everything seems to be good :)"
	echo ""
fi

rm -f $tmp1 $tmp2

echo "Done"

