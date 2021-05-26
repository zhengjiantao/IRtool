# IRtool
A package containing complete intron retention (IR) prediction and quantification

## Requirements
* python >= 3.6
* gensim >= 3.8.3
* pandas >= 1.0.5
* numpy >= 1.18.5
* tensorflow >= 2.0
* samtools >= 1.7
* bedtools = 2.29.2
* iREAD >= 0.8.6
* IRFinder >= 1.2.6
* RSEM >= 1.3.3

## 1. IR prediction
```
# method: iREAD, IRFinder, DeepRetention
bash src/IR_prediction.sh bam_file genome_fasa genome_gtf method
```

## 2. IR quantification
```
# method: RSEM
bash src/IR_quantification.sh bam_file genome_gtf method
```

## 3. Downstream analysis example
```
# downstream: Classification of disease subtypes, single cell
# exp: AD, Cluster
bash src/IR_downstream.sh exp
```
