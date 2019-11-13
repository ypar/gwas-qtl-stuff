#!/bin/bash

# ypar

## for detailed description of each step see enloc github page: https://github.com/xqwen/integrative/tree/master/dev


## run per tissue

rootdir='~/gtex/enloc_v8'
cd ${rootdir}

while read -r tissue; do 

tissuedir=${rootdir}'/dapg/'${tissue}'/'
scriptname='run_enloc_'${tissue}'_bsub'

if [ -s ${scriptname}.sh ]; then
  rm ${scriptname}.sh
fi

cat << EOF >> ${scriptname}.sh
#!/bin/bash
#BSUB -e ${scriptname}.e
#BSUB -o ${scriptname}.o
#BSUB -J ${scriptname}
#BSUB -t 25:00:00
#BSUB -M 50000

cd ${rootdir}

if [ -s ${tissuedir}/eqtl.pip ]; then
  rm ${tissuedir}/eqtl.pip
fi

perl ${rootdir}/integrative/dev/utility/get_max_pip.pl ${tissuedir}/fm_rst > ${tissuedir}/eqtl.pip

if [ -s ${tissuedir}/sig.summary ]; then
  rm ${tissuedir}/sig.summary
fi

perl ${rootdir}/integrative/dev/utility/summarize_dapg_sig.pl ${tissuedir}/fm_rst > ${tissuedir}/sig.summary

if [ -s ${tissuedir}/eqtl.avg.pip ]; then
  rm ${tissuedir}/eqtl.avg.pip
fi

awk '{ total += \$2 } END { print total/NR }' ${tissuedir}/eqtl.pip > ${tissuedir}/eqtl.avg.pip

EOF

done < tissue_list.txt


