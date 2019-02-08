#!/bin/sh

# check sum for fastq files


for i in *.fastq.gz

do 
	md5sum $i >> md5sum.txt

done


