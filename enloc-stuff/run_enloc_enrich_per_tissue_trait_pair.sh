#!/bin/bash

# ypar

## run per tissue trait pair

rootdir='~/gtex/enloc_v8'
cd ${rootdir}

while read -r tissue; do 

  tissuedir=${rootdir}'/dapg/'${tissue}

  while read -r trait; do

    tissuetraitdir=${rootdir}'/enloc/'${tissue}'_'${trait}
    tissuetraitparams=${tissuetraitdir}'/enloc.params'
    mkdir -p ${tissuetraitdir}

    if [ -s ${tissuetraitparams} ]; then
      rm ${tissuetraitparams}
    fi

    cat << EOF >> ${tissuetraitparams}
bin_dir ~/gtex/enloc_v8/bin
gwas_data ~/gtex/v8/gwas/sliced_${trait}_annotated_enloc.txt.gz
qtl_fm_dir ${tissuedir}/fm_rst/
qtl_pip_file ${tissuedir}/eqtl.pip
qtl_summary_file ${tissuedir}/sig.summary
out_dir ${tissuetraitdir}/output
trait ${tissue}_${trait}
EOF

    scriptname=${tissue}_${trait}_enloc_enrich_bsub.sh

    cat << EOF >> ${scriptname}
#!/bin/bash
#BSUB -e ${tissue}_${trait}_enloc_enrich.e
#BSUB -o ${tissue}_${trait}_enloc_enrich.o
#BSUB -J ${tissue}_${trait}_enloc_enrich
#BSUB -M 120000
#BSUB -t 25:00:00

${rootdir}/bin/enloc_enrich ${tissuetraitparams}
chmod -R 755 ${tissuetraitdir}
${tissuetraitdir}/output/prep_mi.cmd
${tissuetraitdir}/output/mi.cmd
${rootdir}/bin/enloc_rcp ${tissuetraitparams}

EOF

    bsub < ${scriptname}
    mv -i ${scriptname} low_level_script/
    sleep 1

  done < trait_list.txt

done < tissue_list.txt


