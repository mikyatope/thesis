
for chr in 2L 2R 3L 3R X
do
	mkdir ${chr}

	for outgroup in  dyak dsim #no_out
	do
		mkdir ${chr}/${outgroup}
		mkdir ${chr}/${outgroup}/7Mb
		mkdir ${chr}/${outgroup}/3Mb
		mkdir ${chr}/${outgroup}/1Mb
		mkdir ${chr}/${outgroup}/500Kb
		mkdir ${chr}/${outgroup}/250Kb

		#N/3
		msa_split ../${chr}.all.${outgroup}.nodupl.fasta -w 7000000,3500000 -r ${chr}.no_out
		mv *.fa ${chr}/${outgroup}/7Mb/
		cat "files moved. Converting to phylip..."
		sh doPhylip.sh ${chr}/${outgroup}/7Mb/
		cat "removing fasta files..."
		rm ${chr}/${outgroup}/7Mb/*.fa


		#N/7
		msa_split ../${chr}.all.${outgroup}.nodupl.fasta -w 3000000,1500000 -r ${chr}.no_out
		mv *.fa ${chr}/${outgroup}/3Mb/
		cat "files moved. Converting to phylip..."
		sh doPhylip.sh ${chr}/${outgroup}/3Mb/
		cat "removing fasta files..."
		rm ${chr}/${outgroup}/3Mb/*.fa

		#N/21
		msa_split ../${chr}.all.${outgroup}.nodupl.fasta -w 1000000,500000 -r ${chr}.no_out
		mv *.fa ${chr}/${outgroup}/1Mb/
		cat "files moved. Converting to phylip..."
		sh doPhylip.sh ${chr}/${outgroup}/1Mb/
		cat "removing fasta files..."
		rm ${chr}/${outgroup}/1Mb/*.fa

		#N/41
		msa_split ../${chr}.all.${outgroup}.nodupl.fasta -w 500000,250000 -r ${chr}.no_out
		mv *.fa ${chr}/${outgroup}/500Kb/
		cat "files moved. Converting to phylip..."
		sh doPhylip.sh ${chr}/${outgroup}/500Kb/
		cat "removing fasta files..."
		rm ${chr}/${outgroup}/500Kb/*.fa

		#N/81
		msa_split ../${chr}.all.${outgroup}.nodupl.fasta -w 250000,125000 -r ${chr}.no_out
		mv *.fa ${chr}/${outgroup}/250Kb/
		cat "files moved. Converting to phylip..."
		sh doPhylip.sh ${chr}/${outgroup}/250Kb/
		cat "removing fasta files..."
		rm ${chr}/${outgroup}/250Kb/*.fa
	done
done
