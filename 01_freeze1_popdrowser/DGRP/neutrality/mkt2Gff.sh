#!/bin/sh

awk '!/^#|^$/ {$3=sprintf("%.4f",$3); print "2L\tmkt_ni\tremark\t"$1"\t"$2"\t"$3"\t.\t.\tName=dgrp_mkt_ni_1kb"}' 2L_mkt.txt > gffs/mkt_ni.gff3 &&
awk '!/^#|^$/ {$4=sprintf("%.4f",$4); print "2L\tmkt_pval\tremark\t"$1"\t"$2"\t"$4"\t.\t.\tName=dgrp_mkt_pval_1kb"}' 2L_mkt.txt > gffs/mkt_pval.gff3 &&

awk '!/^#|^$/ {$3=sprintf("%.4f",$3); print "2R\tmkt_ni\tremark\t"$1"\t"$2"\t"$3"\t.\t.\tName=dgrp_mkt_ni_1kb"}' 2R_mkt.txt >> gffs/mkt_ni.gff3 &&
awk '!/^#|^$/ {$4=sprintf("%.4f",$4); print "2R\tmkt_pval\tremark\t"$1"\t"$2"\t"$4"\t.\t.\tName=dgrp_mkt_pval_1kb"}' 2R_mkt.txt >> gffs/mkt_pval.gff3 &&

awk '!/^#|^$/ {$3=sprintf("%.4f",$3); print "3L\tmkt_ni\tremark\t"$1"\t"$2"\t"$3"\t.\t.\tName=dgrp_mkt_ni_1kb"}' 3L_mkt.txt >> gffs/mkt_ni.gff3 &&
awk '!/^#|^$/ {$4=sprintf("%.4f",$4); print "3L\tmkt_pval\tremark\t"$1"\t"$2"\t"$4"\t.\t.\tName=dgrp_mkt_pval_1kb"}' 3L_mkt.txt >> gffs/mkt_pval.gff3 &&

awk '!/^#|^$/ {$3=sprintf("%.4f",$3); print "3R\tmkt_ni\tremark\t"$1"\t"$2"\t"$3"\t.\t.\tName=dgrp_mkt_ni_1kb"}' 3R_mkt.txt >> gffs/mkt_ni.gff3 &&
awk '!/^#|^$/ {$4=sprintf("%.4f",$4); print "3R\tmkt_pval\tremark\t"$1"\t"$2"\t"$4"\t.\t.\tName=dgrp_mkt_pval_1kb"}' 3R_mkt.txt >> gffs/mkt_pval.gff3 &&

awk '!/^#|^$/ {$3=sprintf("%.4f",$3); print "X\tmkt_ni\tremark\t"$1"\t"$2"\t"$3"\t.\t.\tName=dgrp_mkt_ni_1kb"}' X_mkt.txt >> gffs/mkt_ni.gff3 &&
awk '!/^#|^$/ {$4=sprintf("%.4f",$4); print "X\tmkt_pval\tremark\t"$1"\t"$2"\t"$4"\t.\t.\tName=dgrp_mkt_pval_1kb"}' X_mkt.txt >> gffs/mkt_pval.gff3 &&

#change nan/null values
sed -i 's/NULL/0.00/g' gffs/*.gff3

exit 0
