### Guidelines for CFG Core Personnel: 
    Protocol Guide: Demultiplex NextSeq 500 ---- for mad7 scripts only

    Author: FA
    Version: 1.2

    Last Updated: 04/26/18


    Login to BiX account on the lab's supercomputer

    Browse to Data-hd drive;

    Create a new folder using the following format: (Use the PI's Name)
    FirstName.LastName.dateOfrun (e.g. for a run performed for Dr. Sammons on november 1st 2017, the folder name will be named as Morgan.Sammons.171101)

    Open this folder and copy the Illumina's Run Folder and the SampleSheet.csv file in this directory. 

    Copy and Paste the bcl2fastq.usage.v.1.2.sh file into this Folder as well (The script can be found in BIX flashdrive /crc.scripts/ or on Desktop).

    Open the SampleSheet.csv template (also found in BIX flashdrive /crc.scripts) file and verify the Sample Names and the corresponding indices in the index and index2 columns.
    Make sure to keep the Sample Projects column empty in the sample sheet.  
    SaveAs the SampleSheet.csv into the same directory as the runfolder. (Note: Do not rename the Sample Sheet, leave it as SampleSheet.csv; 
    Also it has to be in the same directory as the run folder, not inside the run folder :) ) 

    Open terminal:
    type cd and hit return
    type cd again and drag the Illumina's Run folder into the terminal window and hit return; this will copy the path of the run folder into the terminal line. 

    type cd .. and hit return

    type ls and hit return; verify that on the terminal window user is able to see the RunFolder Name, SampleSheet.csv and bcl2fastq.usage.sh 

    Once verified:

    type bash bcl2fastq.usage.v.1.2.sh and hit return.

    user will be asked to enter the name of the run folder followed by the PI's name by the terminal

    Enter the name of the run folder (copy and paste just the name) and hit return; Enter the name of the PI in the format specified in the terminal example; 
    The script will then start running and will continue to run until the files are 

    1. demultiplexed and converted to .fastq.gz format; 
    2. fastqc has run and generated its output files;
    3. multiqc output html report has been generated; and 
    4. the fastqs are tar zipped for upload onto server; 

    Once the script is finished the user is able to find all the files in their corresponding folders in the same directory as the run folder. 






