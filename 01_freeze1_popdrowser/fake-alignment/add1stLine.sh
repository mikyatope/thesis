#!/bin/sh
for i in `ls *.raw`
do
        sed 1i\>$i $i > fasta/$i.fasta 
done
exit 0
