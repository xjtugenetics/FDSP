#! /bin/sh

################FILES##################

#all_SNPs="/home/ys/bone_new_ld098/2_get_usresult/all_chr_snp_afeur_001_chr6"
all_SNPs=`grep "all_SNPs_dir" ../prep/files.txt | cut -f2`
positive_set="../data/positive_set.bed"

echo "Loading files"

################PARAMETERs#############

echo "Reading parameters"
rm -f ../data/negative_set_${dist}.tmp


num=`grep "na_distant" ../prep/parameters.txt | cut -f2 | awk -F',' '{print NF-1}'`

echo "Start selecting...this step may takes several hours. Please wait until this job finishes."

################Start#################
for  n in `seq 0 ${num}`;do

#Prepare distant para

list=(`grep "na_distant" ../prep/parameters.txt | cut -f2 | awk -F',' '{for(i=1;i<=NF;i++)print $i}'`)	
dist=${list[${n}]}
#cat ${all_SNPs}/* > ../data/all_SNPs.maf

#File preparation
rm -f ../data/positive.csv.tmp

cat ${positive_set} | while read line; do
SNP=`echo ${line} | awk '{print $4}'`
grep -w "${SNP}" ${all_SNPs}/* | sort | uniq -c | awk '{print $2}'  >> ../data/positive.csv.tmp
done
tmp_positive="../data/positive.csv.tmp"

#Produce a single negative set file
#cat ${positive_set} | while read line; do
cat ${tmp_positive} | while read line; do
#SNP=`echo ${line} | awk '{print $4}'`
#position=`echo ${line} | awk '{print $3}'`
position=`echo ${line} | cut -d ',' -f2`
MAF=`echo ${line} |cut -d ',' -f4`
chr=`echo ${line} |cut -d ',' -f1`
#MAF=`grep -w "${SNP}" ${all_SNPs}/* | cut -d ',' -f4`
cat ${all_SNPs}/* | awk -F',' -v var=${chr} '{if($1==var) print $0}'  | awk -F',' -v var1=${position} -v var2=${MAF} -v var3=${dist} 'function abs(x) {return x <0 ? -x : x} OFS=","{if(abs($2-var1) <= var3 && abs($4-var2) <= 0.05 )print $0}'  >> ../data/negative_set_${dist}.tmp
done


cat ../data/negative_set_${dist}.tmp | sort | uniq -c | awk '{print $2}' > ../data/negative_set_${dist}.list.tmp
f_list[${n}]=../data/negative_set_${dist}.list.tmp
echo "Finished file ${dist}"
done

#############Intersec##################

cat ../ld/*.ld | grep -w -v 'CHR_B' | awk 'OFS="\t"{if($7>0.1)print $4,$5-1,$5,$6}' | sort | uniq -c | awk 'OFS="\t"{print $2,$3,$4,$5}'  > ../ld/ld_0.1.tmp
length=${#list[@]}
for ((i=1; i<=$length; i++))
do
	cat ${f_list[@]:0:${i}} | grep -w -v 'CHR_B' | sort | uniq -u | awk -F',' 'OFS="\t"{print $1,$2-1,$2,$3}' >  ../data/negative_set_${i}.list.tmp
	cat  ../data/negative_set_${i}.list.tmp ../ld/ld_0.1.tmp |  sort | uniq -d  >  ../data/dup.tmp
	cat ../data/negative_set_${i}.list.tmp  ../data/dup.tmp | sort | uniq -u > ../data/negative_set_${i}.list
	echo Finished intersection ${i}	
	done
	

#####################Proportion################
rm -f ../data/negative_set.bed.tmp

po_size=`wc -l ../data/positive_set.bed |awk '{print $1}'`
na_size_t=$[po_size*20]
n=1
for ((i=1; i<=$length; i++))
do	
	na_size=`wc -l  ../data/negative_set_${i}.list | awk '{print $1}'`
	size=`echo "$na_size * $n" | bc`
	n=`echo "$n - 0.2" | bc`
	if [ ${i} = ${length} ]; then
		na_size_p=`wc -l ../data/negative_set.bed.tmp |awk '{print $1}'`
		if [ ${na_size_p} -lt ${na_size_t} ] ;then
		line_sec=`echo "$na_size_t - $na_size_p" | bc`
		cat ../data/negative_set_${i}.list | shuf -n $line_sec | awk -F',' 'OFS="\t"{print "chr"$1,$2,$3,$4}' >> ../data/negative_set.bed.tmp
		fi
	else
	cat ../data/negative_set_${i}.list | shuf -n $size | awk -F',' 'OFS="\t"{print "chr"$1,$2,$3,$4}' >> ../data/negative_set.bed.tmp
	
	fi
done

cat ../data/negative_set.bed.tmp | sort   > ../data/negative_set.bed

rm -f  ../data/*.tmp
rm -f ../ld/ld_0.1.tmp

echo "Nagetive set selection finished"
