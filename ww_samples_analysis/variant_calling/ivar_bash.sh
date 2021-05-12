#!/bin/bash


foundFiles=($(find $1 -name '*.primertrimmed.rg.sorted.bam' -type f)) 

for file in "${foundFiles[@]}"
do
	cd $(dirname $file)
 
	tmp="${file##*/}"
	prefix="${tmp%.primertrimmed.rg.sorted.bam}"
	echo $file
	echo $prefix
	samtools mpileup -aa -A -d 0 -B -Q 0 --reference "/home/diamant/sars-cov2.fa" $file | ivar variants -p $prefix -t 0.03 -r "/home/diamant/sars-cov2.fa" -g /home/diamant/Sars_cov_2.ASM985889v3.101.gff3
done

   
