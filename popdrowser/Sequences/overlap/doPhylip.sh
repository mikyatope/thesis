#!/bin/sh
folder=$1

for i in `ls ${folder}*.fa`
do
	perl fasta2phylip.pl -in $i -out $i.phy &
done
wait
exit 0
