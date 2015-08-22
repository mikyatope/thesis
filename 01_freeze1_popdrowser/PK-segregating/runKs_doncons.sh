#!/bin/sh
date

for size in  100000 50000 10000 #500000 1000000 
do
	for xfold in 0 2 4
	do
		for chr in 2L 2R 3L 3R X 
		do
			perl neutral2Ks_doncons.pl -seq recoded_seq/${chr}_recoded.txt -aln /home/share/${chr}.dsim.dgrp_consensus.fasta -lab ks_4fold_${size} -winsize ${size} -chr ${chr} -fold ${xfold} > results/Ks/${chr}.Ks${xfold}fold.${size}.doncons.gff3 & 
		done
		wait
	done
done

date

exit 0
