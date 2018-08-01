## bixHPC

### Required Packages:
    bcl2fastq
    fastqc
    multiqc
    star
    homer
    DESeq2

## Scripts: 

    barcodeSplit.pl: 
        Demultiplexes based on barcodes of interest.
    
    bcl2fastq.HPC.sh:
            Demultiplexes Illumina NextSeq 500 samples (bcl intensities to fastq).
            Performs quality control assessments using FASTQC. 
            Generates a single html report plotting all quality control metrics.
        
        Usage: 
            $ sbatch bcl2fastq.HPC.sh NameOfRunFolder NameOfSampleSheet NameOfPI         
                           
    star.align.HPC.sh: 
            Aligns PE or SE fastq samples to reference genome using star ultrafast aligner. 
            Generates bedgraph files for UCSC genome browser visualization. 
        
        Usage: 
            $ sbatch star.align.HPC.sh <SE or PE> <ORG> <ASSEMBLY> 
    
    DESeq2.R:
        Performs differential gene expression analysis of RNA-seq samples. 
        (rawCounts used are either from STAR --quantMode GeneCounts or HTSEQ outputs)
    
    reHash.py:
        Rename utility.
    





    
