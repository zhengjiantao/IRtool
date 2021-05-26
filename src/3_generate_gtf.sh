#!/bin/bash
# 令牌start
[ -e ./fd1 ] || mkfifo ./fd1
exec 3<> ./fd1  
rm -rf ./fd1   
for i in `seq 1 20`;
do
    # echo 每次输出一个换行符,也就是一个令牌
    echo >&3              
done
# 令牌end

main="/root/userfolder/quantify"
method="iREAD"
dataset=("Breast_GSE48213")
for GSE in ${dataset[@]};
do
    echo "processing ${GSE}..."
    cd $main/$GSE/$method
    if [ ! -d 'intronGTF' ];then
        mkdir intronGTF
    fi
    cat $main/$GSE/SRR.list | while read line
    do
		read -u3
		{
			if [ $GSE = "Cell_GSE57249" ];
			then
				gtf=$main/Mus_musculus.GRCm38.75.gtf
			else
				gtf=$main/Homo_sapiens.GRCh38.77.gtf
			fi
	        python $main/code/generate_intron_transcript_gtf.py $gtf $line intronGTF
			echo >&3
		}&
    done
	wait
done

exec 3<&-                       # 关闭文件描述符的读
exec 3>&-                       # 关闭文件描述符的写
