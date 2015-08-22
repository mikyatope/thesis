#!/bin/sh

date

#create results directory 
mkdir results/

for j in `ls /home/share/dgrp.freeze1/*.all.no_out.nodupl.fasta.phy`
do
	i=`basename "$j" .all.no_out.nodupl.fasta.phy`
	
	for size in 100000 50000 10000 1000 500 100 50
	do
		#Variation NO-outgroup
		#guardar temps -> 	/usr/bin/time -o logs/$i.$size.$estimate.time -f "%E real \n%U user \n%S sys" 
		for estimate in pi s eta eta_e theta
		do
			perl VariScan2.pl -f $j -o none -RefSeq 2 -chr $i -trmod ${estimate} -w $size -j $size -m 1 -wt 3 -LDsin 0 -CompleteDeletion 0 -FixNum 1 -NumNuc 140 -$estimate -seed 5 &	
		done
		#wait
				
		#Neutrality NO outgroup
		for estimate in tajD fuliDast fuliFast 
		do
			perl VariScan2.pl -f $j -o none -RefSeq 2 -chr $i -trmod ${estimate} -w $size -j $size -m 1 -wt 3 -LDsin 0 -CompleteDeletion 0 -FixNum 1 -NumNuc 140 -$estimate -seed 5  &	
		done
		#wait
		
		#LD NO outgroup
		for estimate in LD_sites D Dabs Dprime Dprimeabs r2 h hd FuFs
		do
			perl VariScan2.pl -f $j -o none -RefSeq 2 -chr $i -trmod ${estimate} -w $size -j $size -m 1 -wt 3 -LDsin 1 -CompleteDeletion 1 -FixNum 1 -NumNuc 140 -$estimate -seed 5  &	
		done
		wait
		
	done
done


for j_out in `ls /home/share/dgrp.freeze1/*.all.dsim.nodupl.fasta.phy`
do
	i_out=`basename "$j_out" .all.dsim.nodupl.fasta.phy`
	
	for size in 100000 50000 10000 1000 500 100 50
	do
		#Variation outgroup
		for estimate in k s_inter
		do
			perl VariScan2.pl -f $j_out -o first -RefSeq 2 -chr $i_out -trmod ${estimate}_dsim -w $size -j $size -m 1 -wt 3 -LDsin 0 -CompleteDeletion 0 -FixNum 1 -NumNuc 140 -$estimate -seed 5  &	
		done

		#Neutrality outgroup
		for estimate in fuliD fuliF FayWuH
		do
			perl VariScan2.pl -f $j_out -o first -RefSeq 2 -chr $i_out -trmod ${estimate}_dsim -w $size -j $size -m 1 -wt 3 -LDsin 0 -CompleteDeletion 0 -FixNum 1 -NumNuc 140 -$estimate -seed 5  &	
		done
		
		#LD with outgroup
		for estimate in LD_sites D Dabs Dprime Dprimeabs r2 h hd FuFs
		do
			perl VariScan2.pl -f $j_out -o first -RefSeq 2 -chr $i_out -trmod ${estimate}_dsim -w $size -j $size -m 1 -wt 3 -LDsin 1 -CompleteDeletion 1 -FixNum 1 -NumNuc 140 -$estimate -seed 5  &	
		done
		wait
		
	done
done


for j_out in `ls /home/share/dgrp.freeze1/*.all.dyak.nodupl.fasta.phy`
do
	i_out=`basename "$j_out" .all.dyak.nodupl.fasta.phy`
	
	for size in 100000 50000 10000 1000 500 100 50
	do
		#Variation outgroup
		for estimate in k s_inter
		do
			perl VariScan2.pl -f $j_out -o first -RefSeq 2 -chr $i_out -trmod ${estimate}_dyak -w $size -j $size -m 1 -wt 3 -LDsin 0 -CompleteDeletion 0 -FixNum 1 -NumNuc 140 -$estimate -seed 5  &	
		done

		#Neutrality outgroup
		for estimate in fuliD fuliF FayWuH
		do
			perl VariScan2.pl -f $j_out -o first -RefSeq 2 -chr $i_out -trmod ${estimate}_dyak -w $size -j $size -m 1 -wt 3 -LDsin 0 -CompleteDeletion 0 -FixNum 1 -NumNuc 140 -$estimate -seed 5  &	
		done
		
		#LD with outgroup
		for estimate in LD_sites D Dabs Dprime Dprimeabs r2 h hd FuFs
		do
			perl VariScan2.pl -f $j_out -o first -RefSeq 2 -chr $i_out -trmod ${estimate}_dyak -w $size -j $size -m 1 -wt 3 -LDsin 1 -CompleteDeletion 1 -FixNum 1 -NumNuc 140 -$estimate -seed 5  &	
		done
		wait
		
	done
done




date

exit 0
