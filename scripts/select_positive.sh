#! /bin/sh

############PARAMETERS##################

po_threshold=`grep "po_threshold" ../prep/parameters.txt | cut -f2`

cat ../ld/*.ld | grep -w -v 'CHR_B' | awk -v var=${po_threshold} '{if($7>var)print "chr"$4,$5-1,$5,$6}' |  sort | uniq -c | awk 'OFS="\t"{print $2,$3,$4,$5}'  > ../data/positive_set.bed



echo 'Finished positive set selection.'
