#!/bin/bash
#
# noise_00761_0000385_0.05_g_2Mpc.rebin.fits
# noise_00761_0000174_0.05_g_2Mpc.rebin.fits
# noise_00761_0000126_0.05_g_2Mpc.rebin.fits
# noise_00761_0000120_0.05_g_2Mpc.rebin.fits
# noise_00761_0000078_0.05_g_2Mpc.rebin.fits
#
path_simulations="/home/ellien/LSST_ICL/simulations/out1/*rebin*"


for file in $( ls $path_simulations )
do
      echo "Launch Dawis on file ${file}"
      qsub qsub_dawis.sh -v ncl=${a:12:7},nf=${file}
      sleep 2
done
