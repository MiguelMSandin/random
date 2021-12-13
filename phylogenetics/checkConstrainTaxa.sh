#!/bin/bash

usage()
{
    echo ""
    echo "Usage: checkConstrainTaxa.sh -d DIRECTORY -o OUTPUT"
    echo ""
    echo "  -d    a directory containing the selected taxa for each taxon."
    echo "  -o    an output file containing the unique summarized taxons at a higher taxonomic level. If exists, will be overwritten."
    echo ""
    echo "Please check the code according "
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

i=0
LENGTH=$(ls -1 $DIRECTORY | wc -l)
for FILE in $(ls $DIRECTORY)
do
	((i=i+1))
	echo -n -e "\r  $i/$LENGTH: $FILE                    "
	
	if [ -s "$DIRECTORY/$FILE" ]; then
		echo "$FILE" >> $OUT
		
		echo "	PR2 search" >> $OUT
		if [ $(grep -c "PR2_" "$DIRECTORY/$FILE") == 0 ]; then
			echo "No matches" >> $OUT
		fi
		grep "PR2_" "$DIRECTORY/$FILE" | \
			cut -d'|' -f 3,4,5,6 | \
			sort | \
			uniq | \
			sed 's/|/\t/g' >> $OUT
			
		echo "	PacBio search" >> $OUT
		if [ $(grep -c "PacBio_" "$DIRECTORY/$FILE") == 0 ]; then
			echo "No matches" >> $OUT
		fi
		grep "PacBio_" "$DIRECTORY/$FILE" | \
			sed 's/_X/-X/g' | \
			sed 's/.*Eukaryota_//g' | \
			cut -d'_' -f 1,2,3,4,5 | \
			sort | \
			uniq | \
			sed 's/_/\t/g' >> $OUT
			
		echo "" >> $OUT
	else
		rm -f "$DIRECTORY/$FILE"
	fi
done

echo -e "\rDone                                        "
