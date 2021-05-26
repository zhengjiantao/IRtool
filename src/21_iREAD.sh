#!/bin/bash

main="/root/userfolder/quantify"

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

dataset=("Cell_GSE57249")
for GSE in ${dataset[@]};
do
    cd $main/$GSE
    if [ ! -d 'iREAD' ];then
        mkdir iREAD
    fi
    cat SRR.list.again | while read line
    do
		echo $line
        read -u3
        {
            bam=${line}_STAR/${line}Aligned.sortedByCoord.out.bam
            if [ $GSE = "AD_GSE125583" -o $GSE = "Cell_GSE57249" ];
            then
                # single
                map_reads=$(samtools flagstat $bam --threads 15 | grep 'total' | cut -d' ' -f1)
            else
                # pair
                map_reads=$(samtools flagstat $bam --threads 15 | grep 'properly' | awk '{print $1}')
            fi

			if [ $GSE = "Cell_GSE57249" ];
			then
				bed_dir="/root/userfolder/software/iread-0.8.9/meta/intron_mouse_3875.bed"
			else
				bed_dir="/root/userfolder/software/iread-0.8.9/meta/intron_human_3877.bed"
			fi

			iread.py $bam $bed_dir -o iREAD -t $map_reads -k 10 -n 10 -j 1 -f 0.1 -e 0.5
            echo >&3
        }&
    done
    wait
done

exec 3<&-                       # 关闭文件描述符的读
exec 3>&-                       # 关闭文件描述符的写
