#!/bin/sh 

# NOTE THE SEQUENE OF FILE NAMES ARE IMPORTANT, PLEASE KEEP TRACK OF THE ORDER OF FILE NAMES AND CONDITIONS
# THE CODE NEEDS TO BE MANUALLY CHANGED EACH TIME FOR A DESEQ2 RUN


FILE1=ATF3_DMSO_1.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE2=ATF3_DMSO_2.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE3=ATF3_DMSO_3.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE4=ATF3_Nutlin_1.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE5=ATF3_Nutlin_2.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE6=ATF3_Nutlin_3.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE7=p53_DMSO_1.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE8=p53_DMSO_2.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE9=p53_DMSO_3.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE10=p53_Nutlin_1.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE11=p53_Nutlin_2.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE12=p53_Nutlin_3.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE13=WT_DMSO_1.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE14=WT_DMSO_2.fq.gz.ReadsPerGene.out.tab.rawCounts
FILE15=WT_DMSO_3.fastq.gz.ReadsPerGene.out.tab.rawCounts
FILE16=WT_Nutlin_1.fastq.gz.ReadsPerGene.out.tab.rawCounts
FILE17=WT_Nutlin_2.fastq.gz.ReadsPerGene.out.tab.rawCounts
FILE18=WT_Nutlin_3.fastq.gz.ReadsPerGene.out.tab.rawCounts


paste $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6 $FILE7 $FILE8 $FILE9 $FILE10 $FILE11 $FILE12 $FILE13 $FILE14 $FILE15 $FILE16 $FILE17 $FILE18 > countMatrix.txt #add more files if necessary



