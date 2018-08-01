## bixHPC

### Required Packages:
    bcl2fastq
    fastqc
    multiqc
    star
    homer
    DESeq2

## Scripts: 

    barcodeSplit.pl     -> Demultiplexes based on barcodes of interest.
    bcl2fastq.HPC.sh    -> Demultiplexes Illumina NextSeq 500 samples (bcl intensities to fastq).
                           Performs quality assessment using FASTQC.
                           Generates a single html (json-based) report plotting all quality metrics.
    star.align.HPC.sh   -> Aligns PE or SE fastq samples to reference genome using star ultrafast aligner. 
    DESeq2.R            -> Differential gene expression analysis of RNA-seq samples.
    reHash.py           -> Rename utility.
    





