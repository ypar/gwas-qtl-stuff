#!/bin/bash

# ypar

# prepare input for enloc
for f in sliced_*_annotated.txt; do
  n=$(echo ${f} | sed 's/annotated.txt/annot4enloc.txt.gz/g')
  awk -F'\t' -v OFS='\t' '{print $8, $4, $9}' ${f} | gzip - > ${n}
done

