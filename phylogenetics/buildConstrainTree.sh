#!/bin/bash

usage()
{
    echo ""
    echo "Usage: buildConstrainTree.sh -d DIRECTORY -t TREE_IN -o TREE_OUT"
    echo ""
    echo "  -d    A directory containing the selected taxa for each taxon."
    echo "  -t    A newick tree."
    echo "  -o    The exported newick tree."
    echo "  -n    If selected, will remove node names."
    echo ""
}

while getopts "hd:t:o:n" opt; do
        case ${opt} in
                h )
                        usage
                        exit
                        ;;
                d )
                        DIRECTORY=$OPTARG
                        ;;
                t )
                        TREE=$OPTARG
                        ;;
                o )
                        OUT=$OPTARG
                        ;;
                n )
                        NODE="TRUE"
                        ;;
        esac
done

if test -f "$OUT"; then
	echo "$OUT file exists, deleting..."
	rm -f $OUT
fi

cp $TREE $OUT

tmp=$(mktemp --tmpdir=$(pwd))
tmp1=$(mktemp --tmpdir=$(pwd))
tmp2=$(mktemp --tmpdir=$(pwd))

i=0
LENGTH=$(ls -1 $DIRECTORY | wc -l)
for FILE in $(ls $DIRECTORY)
do
	((i=i+1))
	echo -n -e "\r  $i/$LENGTH: $FILE                    "
	
	GROUP="${FILE%.*}"
	
	sed "1s/^//" $DIRECTORY/$FILE > $tmp
	sed -i "s/$/INSERTCOMMA/" $tmp
	sed -i "$ s/INSERTCOMMA//" $tmp
	
	if [[ $NODE == "TRUE" ]]; then
		if (( $(wc -l <$tmp) == 1 )); then
			sed "s/$GROUP.*//" $OUT > $tmp1
			sed "s/.*$GROUP//" $OUT > $tmp2
		else
			sed "s/$GROUP.*/(/" $OUT > $tmp1
			sed "s/.*$GROUP/)/" $OUT > $tmp2
		fi
	else
		sed "s/$GROUP.*/(/" $OUT > $tmp1
		sed "s/.*$GROUP/)$GROUP/" $OUT > $tmp2
	fi
	
	{ cat $tmp1; cat $tmp; cat $tmp2; } > $OUT
	
	sed 's/INSERTCOMMA/,/g' $OUT > $tmp
	tr '\n' ' ' < $tmp > $OUT
	
done

sed -i 's/ //g' $OUT
sed -i 's/;/;\n/g' $OUT

rm -f $tmp $tmp1 $tmp2

echo -e "\rDone                                        "
