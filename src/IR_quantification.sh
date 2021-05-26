#!/bin/bash

bam=$1
gtf=$2
method=$3

main=/home/zhengjt/IRtool
export $main

# construct intron-retaining isoforms
python generate_intron_transcript_gtf.py $gtf $bam .
wait

if [ $method = "RSEM" ];
then
	bash 4_RSEM.sh
elif [ $method = "Salmon" ];
then
	echo "Salmon is wait for joining..."
else
	echo "please choose one isoform quantification method!"
fi
