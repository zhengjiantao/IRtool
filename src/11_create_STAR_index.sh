#!/bin/bash

outdir=".."

if [ ! -d ${outdir}/STAR_index ];then
	mkdir ${outdir}/STAR_index
fi

fa_file="../Mus_musculus.GRCm38.75.fa" #"../Homo_sapiens.GRCh38.dna.SORTED.fa"
gtf_file="../Mus_musculus.GRCm38.75.gtf" #"../Homo_sapiens.GRCh38.77.gtf"

STAR \
--runThreadN 30 \
--runMode genomeGenerate \
--genomeDir ${outdir}/STAR_index \
--genomeFastaFiles $fa_file \
--sjdbGTFfile $gtf_file \
--sjdbOverhang 100
