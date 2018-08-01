## bixHPC

### Required Packages:
    bcl2fastq
    fastqc
    multiqc
    star
    homer
    DESeq2

## Scripts: 

    barcodeSplit.pl     -> Demultiplexes based on barcodes of interest
    bcl2fastq.HPC.sh    -> Demultiplexes Illumina NextSeq 500 samples, followed by FASTQC run, generation of multiqc html report and tar zip of samples
    star.align.HPC.sh   -> Alignment to reference genome using star ultrafast aligner
    DESeq2.R            -> Differential gene expression analysis of RNA-seq samples
    reHash.py           -> Rename utility
    





