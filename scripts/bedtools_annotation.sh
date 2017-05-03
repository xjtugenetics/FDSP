#! /bin/sh

################FILES##################

positive_set=../data/positive_set.bed
negative_set=../data/negative_set.bed
annotation=`sed '1d' ../prep/annotation.txt | cut -f2`
annotation_file=`grep "annotation_file" ../prep/files.txt | cut -f2`
title=`sed '1d' ../prep/annotation.txt | cut -f1`

################Annotation################
rm -f ../data/positive_set_annotation.bed
rm -f  ../data/positive_set_annotation.bed

if [ -f "$positive_set" ]; then 
echo  ${title} | sed 's/ /,/g' | awk 'OFS=","{print "SNP",$0,"status"}' > ../data/positive_set_annotation.csv
bedtools annotate -i ${positive_set} -files ${annotation}  -counts | awk  '{for(i=1;i<=3;i++)$i=""; print $0}' | awk  'OFS=","{$(NF+1)="1";print $0}' >> ../data/positive_set_annotation.csv
fi

if [ -f "$negative_set" ]; then
echo  ${title} | sed 's/ /,/g' | awk 'OFS=","{print "SNP",$0,"status"}' > ../data/negative_set_annotation.csv
 bedtools annotate -i ${negative_set} -files ${annotation}  -counts | awk  '{for(i=1;i<=3;i++)$i=""; print $0}' | awk  'OFS=","{$(NF+1)="1";print $0}' >> ../data/negative_set_annotation.csv
 fi

if [ ! $annotation_file = "NA" ]; then
	echo  ${title} | sed 's/ /,/g' | awk 'OFS=","{print "SNP",$0}' > ../data/annotation_file.csv
	 bedtools annotate -i ${annotation_file} -files ${annotation}  -counts | awk  '{for(i=1;i<=3;i++)$i=""; print $0}' |   awk  'OFS=","{print $0}' >> ../data/annotation_file.csv
 fi
