#!/bin/sh
date



for fixnum in 140 130 120 110
do
	perl neutral3Ps.pl -seq recoded_seq/2L_recoded.txt -aln /home/share/2L.all.no_out.nodupl.fasta -lab ps_4fold_100000_fixnum${fixnum} -winsize 100000 -chr 2L -fold 4 -fixnum ${fixnum} > results/2L.Ps.4fold.100000.fixnum${fixnum}.countdown.gff3 & 
done
wait

for fixnum in 145 135 125 115 105
do
	perl neutral3Ps.pl -seq recoded_seq/2L_recoded.txt -aln /home/share/2L.all.no_out.nodupl.fasta -lab ps_4fold_100000_fixnum${fixnum} -winsize 100000 -chr 2L -fold 4 -fixnum ${fixnum} > results/2L.Ps.4fold.100000.fixnum${fixnum}.countdown.gff3 & 
done
wait

date

exit 0
