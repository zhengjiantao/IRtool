#!/bin/bash

bam=$1
genome=$2
gtf=$3
method=$4

main=/home/zhengjt/IRtool
export $main

if [ $method = "iREAD" ];
then
	bash 21_iREAD.sh $bam $genome $gtf
elif [ $method = "IRFinder" ];
then
	bash 23_IRFinder.sh $bam $genome $gtf
elif [ $method = "DeepRetention" ];
then
	python 22_DeepRetention.py $bam $genome $gtf
else
	echo "please choose one IR detection method!"
fi
