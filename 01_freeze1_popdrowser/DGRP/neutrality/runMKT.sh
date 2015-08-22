/usr/bin/time -o mkt-2L_10k.time -f "%E real \n%U user \n%S sys" perl neutral.pl -aln ../datasets/with_dsim/x4/2L_w_dsim_x4.fasta  -seq ../recode_seq/2L_recoded.txt  -winsize 10000 > 2L_10k.mkt &
/usr/bin/time -o mkt-2R_10k.time -f "%E real \n%U user \n%S sys" perl neutral.pl -aln ../datasets/with_dsim/x4/2R_w_dsim_x4.fasta  -seq ../recode_seq/2R_recoded.txt  -winsize 10000 > 2R_10k.mkt &
/usr/bin/time -o mkt-3L_10k.time -f "%E real \n%U user \n%S sys" perl neutral.pl -aln ../datasets/with_dsim/x4/3L_w_dsim_x4.fasta  -seq ../recode_seq/3L_recoded.txt  -winsize 10000 > 3L_10k.mkt &
/usr/bin/time -o mkt-3R_10k.time -f "%E real \n%U user \n%S sys" perl neutral.pl -aln ../datasets/with_dsim/x4/3R_w_dsim_x4.fasta  -seq ../recode_seq/3R_recoded.txt  -winsize 10000 > 3R_10k.mkt &
/usr/bin/time -o mkt-X_10k.time -f "%E real \n%U user \n%S sys" perl neutral.pl -aln ../datasets/with_dsim/x4/X_w_dsim_x4.fasta  -seq ../recode_seq/X_recoded.txt  -winsize 10000 > X_10k.mkt &
