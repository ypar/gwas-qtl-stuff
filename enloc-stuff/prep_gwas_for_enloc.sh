#!/bin/bash

# ypar

# run per trait downloaded from gcp

for f in imputed_*.txt.gz;
  do
  trait=$(echo ${f} | sed 's/imputed_//g; s/.txt.gz//g')
  slicefile='sliced_'${trait}'.txt'
  annotated='sliced_'${trait}'_annotated.txt'

  cat << EOF >> prep_gwas_for_enloc_${trait}.sh
#!/bin/bash
#BSUB -e prep_gwas_for_enloc_${trait}.e
#BSUB -o prep_gwas_for_enloc_${trait}.o
#BSUB -J prep_gwas_for_enloc_${trait}
#BSUB -t 25:00:00
#BSUB -M 50000

cd ~/gtex/v8/gwas/
# minus 1 for null position - bedtools bed is 0-based
if [ ! -s ${slicefile} ]; then
  zless -s ${f} | awk -F'\t' -v OFS='\t' '(NR>1){print \$3,\$4-1,\$4,\$2,\$10}' > ${slicefile}
fi
if [ ! -s ${annotated} ]; then
  bedtools intersect -wa -wb -a ldregion.txt -b ${slicefile} > ${annotated}
fi

EOF

done



