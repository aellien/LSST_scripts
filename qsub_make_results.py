#!/bin/bash
#PBS -o /home/ellien/LSST_ICL/logs/make_results.out
#PBS -j oe
#PBS -N maker
#PBS -l nodes=1:ppn=3,walltime=47:00:00
#PSB -S /bin/bash

module load intelpython/3-2020.4
python /home/ellien/LSST_ICL/LSST_scripts/make_results.py


exit 0
