for outgroup in dsim dyak
do
	for chr in 2L 2R 3L 3R X
	do 
		sed -n 7,8p ../${chr}.alloutgroups.nodupl.fasta > ${chr}.${outgroup}.fasta 
		cat ${chr}.${outgroup}.fasta ${chr}.dgrp_consensus.fasta > ${chr}.${outgroup}.dgrp_consensus.fasta 
	done
done
