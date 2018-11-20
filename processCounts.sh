#!/bin/sh

if [ "$1" = "help" ]
    then
    echo ""
    echo " usage = bash processCounts.sh <1> or <2> ; where 1 is for first strand counts and 2 is for reverse strand counts "
    echo ""


elif [ "$1" = "1" ]
    then
        for i in *.ReadsPerGene.out.tab
        do
        awk 'NR > 4 {print $1 "\t" $3}' $i > $i.rawCounts
        done
elif [ "$1" = "2" ]
    then
        for i in *.ReadsPerGene.out.tab
        do
        awk 'NR > 4 {print $1 "\t" $4}' $i > $i.rawCounts
        done
else
    for i in *.ReadsPerGene.out.tab
    do
    awk 'NR > 4 {print $1 "\t" $2}' $i > $i.rawCounts
    done
fi












