#!/bin/bash


foundFiles=($(find $1 -name '*.primertrimmed.rg.sorted.bam' -type f)) 

for file in "${foundFiles[@]}"
do
	cd $(dirname $file)
 
	tmp="${file##*/}"
	prefix="${tmp%.primertrimmed.rg.sorted.bam}"
	echo $file
	echo $prefix
	bedtools coverage -abam $file -b "/home/diamant/pool1.bed" > $prefix.depth1.txt
	bedtools coverage -abam $file -b "/home/diamant/pool2.bed" > $prefix.depth2.txt	
done

   
