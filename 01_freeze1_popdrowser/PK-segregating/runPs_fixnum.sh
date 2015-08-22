#!/bin/sh
date

for size in  100000 50000 10000 #500000 1000000 
do
	for xfold in 4 0 2
	do
			for chr in 2L 2R 3L 3R X
			do
				perl neutral3Ps.pl -seq recoded_seq/${chr}_recoded.txt -aln /home/share/${chr}.all.no_out.nodupl.fasta -lab ps_${xfold}fold_${size}_fixnum140 -winsize ${size} -chr ${chr} -fold ${xfold} -fixnum 140 > results/${chr}.Ps.${xfold}fold.${size}.fixnum140.gff3 & 
			done
			wait
	done
done


date

exit 0
