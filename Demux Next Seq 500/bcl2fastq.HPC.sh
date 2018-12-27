#!/bin/sh -l

#  bcl2fastq.HPC.sh
#  BiX.alpha
#
#  Created by Faraz Ahmed on 5/29/18.
#
# initial commit FA180529

#SBATCH -J bcl2fastq.HPC
#SBATCH -o out.bcl2fastq
#SBATCH -e err.bcl2fastq
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000
# To display help run bash bcl2fastq.HPC.sh help

if [  $1 == "help"  ]; then
        echo ""
        echo "#-----------------------------------------------------------------------------------------------------------------------"
        echo "# Run this script from the same directory as the run folder"
        echo "# To execute the script use the following command : sbatch bcl2fastq.HPC.sh NameOfRunFolder NameOfSampleSheet NameOfPI"
        echo "# If the script gives an execution error try the following command chmod 755 bcl2fastq.HPC.sh and re-run"
        echo "#-----------------------------------------------------------------------------------------------------------------------"
        echo ""
else

#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DEMUX
start=`date +%s`

bcl2fastq \
--runfolder-dir $1 \
--output-dir Output.Fastqs \
--sample-sheet $2 \
-r 2 \
-p 8 \
-w 2 \
--no-lane-splitting

#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PROCESS UNDETERMINED
cd Output.Fastqs
mkdir -p ../Undetermined
mv Undetermined* ../Undetermined
mv Reports ../Undetermined
mv Stats ../Undetermined
cd ..

#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FASTQC
mkdir Fastqc.Results

cd Output.Fastqs
for i in *.fastq.gz
do
fastqc -o ../Fastqc.Results -t 4 "$i"
done
cd ..

#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MULTIQC
multiqc -f -n $3.multiqc.results Fastqc.Results
rm -r $3.multiqc.results_data

#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TAR-zip
tar -zcvf $3.tar.gz Output.Fastqs Fastqc.Results


end=`date +%s`
runtime=$((end-start))
echo 'runtime seconds '"$runtime" > runtime.txt

fi

