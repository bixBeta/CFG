#!/bin/sh -l

#  star.align.HPC.sh
#  BiX.alpha
#
#  Created by Faraz Ahmed on 3/8/18.
#
# initial commit FA180608

#SBATCH -J star.align.HPC
#SBATCH -o out.star.align.HPC
#SBATCH -e err.star.align.HPC
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000


if [  $1 == "help" ]; then

	echo ""
	echo "----------------------------------------------------------------------------------------------------------------------"
	echo "navigate to the folder contatining all fastqs and execute using sbatch star.align.HPC.sh <SE or PE> <ORG> <ASSEMBLY>"
    echo "----------------------------------------------------------------------------------------------------------------------"
    echo ""
    echo " 		Specify SE for Single-end and PE for Paired-end libraries"
    echo " 		Available ORG = Homo_sapiens, Mus_musculus, Drosophila_melanogaster, Rattus_norvegicus"
    echo " 		Available ASSEMBLY = hg19, hg38, mm10, dm3, dm6, rn6 "
    echo ""

else

	start=`date +%s`

	if [ "$1" = "PE" ]
	    then
	        ls -1 *R1* > Read1.list
	        ls -1 *R2* > Read2.list
	        paste -d " " Read1.list Read2.list > Reads.list

	elif [ "$1" = "SE" ]
	    then
	        ls -1 *.gz > Reads.list

	else
	    echo "please input appropriate library (SE or PE)"
	fi


	readarray fastqs < Reads.list # creates an array of fastq samples

	for i in "${fastqs[@]}" 
	do
	    FASTQ="$i"
	    DIR=/network/rit/lab/ahmedlab/genomes/ucsc/$2/$3/$3
	    iSUB=`echo "$i" | cut -d ' ' -f1`

	    STAR \
	    --runThreadN 12 \
	    --genomeDir $DIR \
	    --readFilesIn $FASTQ \
	    --readFilesCommand zcat \
	    --outSAMstrandField intronMotif \
	    --outFilterIntronMotifs RemoveNoncanonical \
	    --outSAMtype BAM SortedByCoordinate \
	    --outFileNamePrefix $iSUB. \
	    --limitBAMsortRAM 61675612266 \
	    --quantMode GeneCounts

	done

	multiqc -f -n STAR.multiqc.report .

	mkdir STAR.COUNTS
	mv *.ReadsPerGene.out.tab STAR.COUNTS

	mkdir BAMS
	mv *.bam BAMS

	mkdir STAR.LOGS
	mv *.out *.tab *_STARtmp *.list STAR.LOGS
	mv STAR.multiqc.report_data STAR.LOGS

	cd BAMS 

	    for b in *.bam
	        do 
	        makeTagDirectory "$b".tag.dir "$b"
	    done

	        mkdir TAGDIRS
	        mv *.tag.dir TAGDIRS

	        cd TAGDIRS 
	        for t in *.tag.dir
	        do 
	        makeUCSCfile "$t" -o auto -fsize 5e7 -res 1 -color 106,42,73 -style rnaseq -strand both
	        done

	        mkdir BEDFILES
	        for b in *.tag.dir
	        do
	            cd $b
	            mv *.bedGraph.gz ../BEDFILES
	            cd ..
	        done


	        mv BEDFILES ../../
	        cd ..

	cd ..

	end=`date +%s`
	runtime=$((end-start))
	echo 'runtime seconds '"$runtime" > runtime.align.txt

fi
