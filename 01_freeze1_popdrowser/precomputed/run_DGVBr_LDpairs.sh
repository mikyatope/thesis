#!/bin/sh

date

#create results directory 
mkdir results/LDpairs

for j in `ls /home/share/dgrp.freeze1/*.all.no_out.nodupl.fasta.phy`
do
	i=`basename "$j" .all.no_out.nodupl.fasta.phy`
	
	for size in 100000 
	do
		#Variation NO-outgroup
		#guardar temps -> 	/usr/bin/time -o logs/$i.$size.$estimate.time -f "%E real \n%U user \n%S sys" 
		for estimate in r2
		do
			perl VariScan2.pl -f $j -o none -RefSeq 2 -chr $i -trmod ${estimate} -w $size -j $size -m 1 -wt 3 -LDsin 1 -CompleteDeletion 1 -FixNum 1 -NumNuc 140 -$estimate -seed 5 -LDF ${i}.LDpairs.csv &	
		done
		wait
	done
done


date

exit 0
