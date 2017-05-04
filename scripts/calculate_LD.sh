#! /bin/sh

#################LOAD FILES###################
lead_SNPs=`grep "lead_SNPs" ../prep/files.txt | cut -f2`
bfile=`grep -w "bfile" ../prep/files.txt | cut -f2`
bfile_list=`grep -w "bfile_list" ../prep/files.txt | cut -f2`
r2=`grep "r2" ../prep/parameters.txt | cut -f2`

###############LD###########################

if [ ${bfile} = "NA" ]; then
	cat ${bfile_list} | while read line; do 
	file=${line}
	file_n=`echo $file |awk -F '/' '{print $NF}'`

	echo 'Calculating LD'
	plink --bfile ${file} --r2  --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 ${r2} --ld-snp-list ${lead_SNPs} --out ../ld/${file_n}
	done
else
	file=${bfile}
	file_n=`echo $file |awk -F '/' '{print $NF}'`
	echo ${file_n}
	plink --bfile ${file} --r2  --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 ${r2} --ld-snp-list ${lead_SNPs} --out ../ld/${file_n}
fi
echo 'Calculation finished'
