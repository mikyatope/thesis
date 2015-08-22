#!/bin/sh
date

for size in  50000 10000 #500000 1000000 100000
do
	for chr in 2L 2R 3L 3R X 
	do
		perl createConsensus.pl -seq recoded_seq/${chr}_recoded.txt -aln /home/share/${chr}.all.dsim.nodupl.fasta -lab ks_4fold_${size} -winsize ${size} -chr ${chr} > ${chr}.dsim.dgrp_consensus.fasta & 
	done
	wait
done

date

exit 0