#!/bin/sh

INPUT=$1

OUTDIR=~/NLM_Reproducibility_Workshop/tb_and_arthritis/working
NAME=$(basename $INPUT .txt.gz)
UNZIP=$(basename $INPUT .gz)
OUTFILE=${OUTDIR}/${NAME}_networkanalyst

gunzip -c $INPUT > ${OUTDIR}/${UNZIP}

CLASS_LINE="#CLASS"; 
for i in $(grep '^[^!]' ${OUTDIR}/${UNZIP} | head -1 | awk '{for (i=2; i<=NF; ++i) print $i}'| sed -e 's/"//g');do 
    CLASS=$(grep $i ~/NLM_Reproducibility_Workshop/tb_and_arthritis/data/samples_mod.txt | awk '{print $2}')
    CLASS_LINE="${CLASS_LINE}\t${CLASS}"
done

grep '^[^!]' ${OUTDIR}/${UNZIP} | head -1 | sed 's/^/#/' > /${OUTFILE}.txt
echo ${CLASS_LINE}>> ${OUTFILE}.txt
echo ${CLASS_LINE}
grep '^[^!]' ${OUTDIR}/${UNZIP} | tail -n +2 >> ${OUTFILE}.txt


awk -F"\t"  'NR==2 { for(i=2;i<=NF;i++){if($i){print i}}}' ${OUTFILE}.txt > ${OUTFILE}-col.txt
awk -F"\t" '{print $1}' ${OUTFILE}.txt > ${OUTFILE}_filter.txt
for i in $(cat ${OUTFILE}-col.txt);do awk -F"\t" -v col=$i '{print $col}' ${OUTFILE}.txt > ${OUTFILE}-${i}.txt; paste ${OUTFILE}_filter.txt ${OUTFILE}-${i}.txt > ${OUTFILE}-tmp.txt;mv ${OUTFILE}-tmp.txt ${OUTFILE}_filter.txt; done

rm ${OUTFILE}-*.txt





