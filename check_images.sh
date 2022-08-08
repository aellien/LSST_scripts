#!/bin/bash
#
# 00761_0000001_0.05_xy_r_4Mpc_mu30.rebin.norm.fits
# 00761_0000001_0.05_yz_r_4Mpc_mu30.rebin.norm.fits
#
path_simulations="/home/ellien/LSST_ICL/simulations/out4/"

for dir in TNG-100/
do
  for file in $path_simulations$dir*.fits
  do
    n=$(basename "$file")
    echo $file
    ds9 ${file} \
    -lock scalelimits yes \
    -lock scale yes \
    -lock frame image \
    -lock colorbar yes \
    -scale mode 99.5
  done
done
