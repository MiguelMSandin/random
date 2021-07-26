#!/bin/bash

for i in $(ls | grep ^RAxML); do
     TYPE=$(echo $i | sed -e "s/\..*//g")
     FILE=$(echo $i | sed -e "s/^[^\.]*\.//g")
     NEW=$FILE"."$TYPE
          if [[ $(echo $NEW | grep "_info") ]]; then
        mv $i $NEW".txt"
        echo $i renamed to: $NEW".txt"
        else
        mv $i $NEW".tre"
        echo $i renamed to: $NEW".tre"
        fi
     if [[ $1 = "-c" ]] || [[ $1 = "-clean" ]]; then
        if [[ $(echo $NEW | grep "_bestTree") ]] || [[ $(echo $NEW | grep "_bipartitionsBranchLabels") ]] || [[ $(echo $NEW | grep "_bootstrap") ]]; then 
            rm $NEW
            echo $NEW deleted
        fi
    fi
done
