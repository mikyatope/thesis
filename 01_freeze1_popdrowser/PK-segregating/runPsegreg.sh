#!/bin/sh
date

for size in  50000 100000 10000 #500000 1000000 
do
	for xfold in 4 0 2
	do
		for chr in 2L 2R 3L 3R X
		do
			perl neutral3Seg.pl -seq recoded_seq/${chr}_recoded.txt -aln /home/share/dgrp.freeze1/${chr}.all.no_out.nodupl.fasta -lab Pseg_${xfold}fold_${size}_fixnum140 -winsize ${size} -chr ${chr} -fold ${xfold} -fixnum 140 -low 1 > results/${chr}.Pseg.${xfold}fold.${size}.fixnum140.csv & 
		done
		wait
	done
done


date

exit 0
