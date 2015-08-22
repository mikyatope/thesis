#!/bin/sh
date

for size in  100000 50000 10000 #500000 1000000 
do
	for xfold in 0 4 2
	do
		for discard in 8 #58
		do
			for chr in 2L 2R 3L 3R X
			do
				perl neutral3Ps.pl -seq recoded_seq/${chr}_recoded.txt -aln /home/share/${chr}.all.no_out.nodupl.fasta -lab ks_${xfold}fold_${size} -winsize ${size} -chr ${chr} -fold ${xfold} -discard ${discard} > results/${chr}.Ps.${xfold}fold.${size}.gff3 & 
			done
			wait
		done
	done
done

date

exit 0
