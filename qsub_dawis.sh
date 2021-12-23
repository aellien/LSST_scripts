#!/bin/bash
#PBS -o /home/ellien/LSST_ICL/logs/icl_LSST_${ncl}.out
#PBS -j oe
#PBS -N icl_LSST
#PBS -l nodes=1:ppn=4,walltime=96:00:00
#PSB -S /bin/bash

module load intelpython/3-2020.4
python /home/ellien/Euclid_ICL/scripts/dawis_LSST.py ${nf}


exit 0
