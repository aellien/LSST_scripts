#!/bin/bash
#PBS -j oe
#PBS -N maker
#PBS -l nodes=1:ppn=1,walltime=24:00:00
#PSB -S /bin/bash

module load intelpython/3-2020.4
python /home/ellien/LSST_ICL/LSST_scripts/make_mocks_fits.py


exit 0
