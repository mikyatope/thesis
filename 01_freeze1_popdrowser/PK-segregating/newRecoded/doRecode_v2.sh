#!/bin/sh
date
perl get4fold_v2.pl -cds /home/dgrp/ref_data/cds_recode/5.13/2L_recode.gff -seq /home/dgrp/ref_data/dmel-all-r5.13/chr/dmel-2L-chromosome-r5.13.fasta > newRecoded/2L_recoded.txt &
perl get4fold_v2.pl -cds /home/dgrp/ref_data/cds_recode/5.13/2R_recode.gff -seq /home/dgrp/ref_data/dmel-all-r5.13/chr/dmel-2R-chromosome-r5.13.fasta > newRecoded/2R_recoded.txt &
perl get4fold_v2.pl -cds /home/dgrp/ref_data/cds_recode/5.13/3L_recode.gff -seq /home/dgrp/ref_data/dmel-all-r5.13/chr/dmel-3L-chromosome-r5.13.fasta > newRecoded/3L_recoded.txt &
perl get4fold_v2.pl -cds /home/dgrp/ref_data/cds_recode/5.13/3R_recode.gff -seq /home/dgrp/ref_data/dmel-all-r5.13/chr/dmel-3R-chromosome-r5.13.fasta > newRecoded/3R_recoded.txt &
perl get4fold_v2.pl -cds /home/dgrp/ref_data/cds_recode/5.13/X_recode.gff -seq /home/dgrp/ref_data/dmel-all-r5.13/chr/dmel-X-chromosome-r5.13.fasta > newRecoded/X_recoded.txt &
date
