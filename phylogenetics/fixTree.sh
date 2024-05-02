#! /bin/bash

FILE=$1

if [ $(grep -c "\];$" $FILE) -ge 1 ]; then
	sed -i -E 's/\];/\];\nEnd;/g' $FILE
else
	sed -i -E 's/\);/\);\nEnd;/g' $FILE
fi

awk -i inplace '!seen[$0]++' $FILE
