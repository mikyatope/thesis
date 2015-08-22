#!/bin/sh
date
perl get4fold.pl -cds ../ref_data/cds_recode/5.13/2L_recode.gff -seq ../ref_data/dmel-2L-chromosome-r5.13.fasta > 5.13/2L_recoded.txt &
perl get4fold.pl -cds ../ref_data/cds_recode/5.13/2R_recode.gff -seq ../ref_data/dmel-2R-chromosome-r5.13.fasta > 5.13/2R_recoded.txt &
perl get4fold.pl -cds ../ref_data/cds_recode/5.13/3L_recode.gff -seq ../ref_data/dmel-3L-chromosome-r5.13.fasta > 5.13/3L_recoded.txt &
perl get4fold.pl -cds ../ref_data/cds_recode/5.13/3R_recode.gff -seq ../ref_data/dmel-3R-chromosome-r5.13.fasta > 5.13/3R_recoded.txt &
perl get4fold.pl -cds ../ref_data/cds_recode/5.13/X_recode.gff -seq ../ref_data/dmel-X-chromosome-r5.13.fasta > 5.13/X_recoded.txt &
date