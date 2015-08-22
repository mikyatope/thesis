for chr in 2L 2R 3L 3R X
do
        perl neutral3SegCountSNPs.pl -aln /home/share/dgrp.freeze1/${chr}.all.no_out.nodupl.fasta -seq newRecoded/${chr}_recoded.txt -chr ${chr} > SNPs/${chr}.snpscall.csv &     
done
