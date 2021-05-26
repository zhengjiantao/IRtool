#!/bin/bash

main="/root/userfolder/quantify"
dataset=("Breast_GSE48213")

# 令牌start
[ -e ./fd1 ] || mkfifo ./fd1
exec 3<> ./fd1  
rm -rf ./fd1   
for i in `seq 1 3`;
do
    # echo 每次输出一个换行符,也就是一个令牌
    echo >&3                   
done
# 令牌end

for GSE in ${dataset[@]};
do
    cd $main/$GSE/iREAD
    if [ ! -d 'RSEM' ];then
        mkdir RSEM
    fi
    cat $main/$GSE/SRR.list | while read line
    do
        gtf="intronGTF/${line}intron.gtf"
        read -u3
        {
            #### 1. build reference ####
            ##output: mouse_ref.transcripts.fa
            if [ ! -d "RSEM/${line}_ref" ];then
                mkdir RSEM/${line}_ref
            fi
			if [ $GSE = "Cell_GSE57249" ];
			then
				fa=$main/Mus_musculus.GRCm38.75.fa
			else
				fa=$main/Homo_sapiens.GRCh38.dna.SORTED.fa
			fi
            rsem-prepare-reference --gtf $gtf \
                                --star \
                                -p 10 \
                                $fa \
                                RSEM/${line}_ref/$line
            #### 2. quantify ####
            if [ ! -d "RSEM/${line}_res" ];then
                mkdir RSEM/${line}_res
            fi
            if [ $GSE = "AD_GSE95587" -o $GSE = "Breast_GSE48213" ];
            then
                # pair end
                rsem-calculate-expression -p 10 --paired-end \
                              --star \
                              --estimate-rspd \
                              --append-names \
                              --output-genome-bam \
                              ../${line}_1.fastq ../${line}_2.fastq \
                              RSEM/${line}_ref/$line RSEM/${line}_res/$line
            else
                # single end
                rsem-calculate-expression -p 10 \
                              --star \
                              --estimate-rspd \
                              --append-names \
                              --output-genome-bam \
                              ../${line}.fastq \
                              RSEM/${line}_ref/$line RSEM/${line}_res/$line
            fi
            echo >&3
        }&
    done
    wait
done

exec 3<&-                       # 关闭文件描述符的读
exec 3>&-                       # 关闭文件描述符的写
