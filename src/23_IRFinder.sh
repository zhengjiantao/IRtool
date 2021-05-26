#!/bin/bash

path=$main/data
cd $path

# build REF for mouse and human
# IRFinder -m BuildRefProcess -r IRFinder_REF/Human3877 &
# IRFinder -m BuildRefProcess -r IRFinder_REF/Mouse3875 &

for GSE in ${dataset[@]};
do
    cd $path/$GSE
    if [ ! -d 'IRFinder' ];then
        mkdir IRFinder
    fi
    
    for bam in $(ls *.sortedByName.bam)
    do
        # bam must be unsort
        if [ $GSE = "human_GM12878" ];
        then
            ref=$path/IRFinder_REF/Human3877
        else
            ref=$path/IRFinder_REF/Mouse3875
        fi
        
        IRFinder -m BAM -r $ref \
                 -d IRFinder/$bam \
                 $bam &
        python ${main}/code/IRFinder_determin_label.py \
               IRFinder/$bam/IRFinder-IR-nondir.txt \
               0.1 \
               IRFinder/${bam%.bam}.ir &
    done
    wait
done