#!/bin/sh

mkdir results/variation results/neutral results/ld results/variation/wig results/neutral/wig results/ld/wig results/ld/outgroup results/ld/outgroup/wig


#Variation
for file in *.pi.* *.s.* *.eta.* *.eta_e.* *.theta.* *.k_* *.s_inter_*
do
	mv results/${file} results/variation/
done
mv results/variation/*.wig results/variation/wig/


#Neutrality
for file in *.tajD.* *.fuliDast.* *.fuliFast.* *.fuliD_* *.fuliF_* *.FayWuH_*
do
	mv results/${file} results/neutral/
done
mv results/neutral/*.wig results/neutral/wig/


#LD
for file in *.LD_sites* *.D* *.r2* *.h* *.hd* *.FuFs*    #al loro tots els D Dabs Dprime Dprimeabs
do
	mv results/${file} results/ld/
done
mv results/ld/*_dsim* results/ld/outgroup/
mv results/ld/*_dyak* results/ld/outgroup/
mv results/ld/*.wig results/ld/wig/
mv results/ld/outgroup/*.wig results/ld/outgroup/wig/


exit 0
