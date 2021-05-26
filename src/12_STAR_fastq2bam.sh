#!/bin/bash

[ -e ./fd1 ] || mkfifo ./fd1
exec 3<> ./fd1  
rm -rf ./fd1   

for i in `seq 1 2`;
do
    # echo 每次输出一个换行符,也就是一个令牌
    echo >&3                   
done

main="/root/userfolder/quantify"
index="STAR_human_index"

GSE="Cell_GSE57249"
cat $main/$GSE/SRR.list | while read line
do
	read -u3
	{
		read=$main/$GSE/${line}.fastq
		output=$main/$GSE/${line}_STAR
		if [ ! -d $output ];then
			mkdir $output
		fi
		STAR \
		--runThreadN 15 \
		--readFilesIn ${read} \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix $output/${line} \
		--sjdbFileChrStartEnd $main/$index/sjdbList.out.tab \
		--genomeDir $main/$index
	
		echo >&3
	}&
done
wait

exec 3<&-                       # 关闭文件描述符的读
exec 3>&-                       # 关闭文件描述符的写
