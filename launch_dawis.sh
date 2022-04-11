#!/bin/bash
#
# HorizonAGN  Hydrangea  Magneticum  TNG-100

path="/n03data/ellien/LSST_ICL/simulations/out2"

for dir in TNG-100
do
    #for file in $path/$dir/*norm.fits
    for file in $path/$dir/00099_0000005_0.05_xz_r_4Mpc_mu30.rebin.norm.fits $path/$dir/00099_0000004_0.05_yz_r_4Mpc_mu30.rebin.norm.fits $path/$dir/00099_0000004_0.05_xy_r_4Mpc_mu30.rebin.norm.fits 
    do
          echo "Launch Dawis on file ${file}"
          n=$(basename "$file")
          qsub qsub_dawis.sh -v ncl=${n:0:21},nf=${n},dir=${dir}
          sleep 2
    done
done
