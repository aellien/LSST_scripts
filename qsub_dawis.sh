#!/bin/bash
#PBS -o /home/ellien/LSST_ICL/logs/icl_LSST_${dir}_${ncl}.out
#PBS -j oe
#PBS -N icl_LSST
#PBS -l nodes=1:ppn=1,walltime=24:00:00
#PSB -S /bin/bash

module load intelpython/3-2020.4
echo ${ncl}
echo ${nf}
#python /home/ellien/LSST_ICL/LSST_scripts/dawis_LSST.py ${nf} ${dir}
python /home/ellien/LSST_ICL/LSST_scripts/dawis_LSST_TNG.py ${nf} ${dir}


exit 0
