#!/bin/bash
#
# 00761_0000001_0.05_xy_r_4Mpc_mu30.rebin.norm.fits
# 00761_0000001_0.05_yz_r_4Mpc_mu30.rebin.norm.fits
#
path_simulations="/home/ellien/LSST_ICL/simulations/out2/"
path_wavelets="/home/ellien/LSST_ICL/wavelets/out2/"
run="run1/"

for dir in Magneticum/ Hydrangea/  # #Hydrangea/
do
  for file in $path_simulations$dir/*norm.fits
  do
    n=$(basename "$file")
    echo $file
    ds9 ${file} \
    ${path_wavelets}${dir}${run}${n:0:-5}.results.gal.hardsepBCG.fits\
    ${path_wavelets}${dir}${run}${n:0:-5}.results.iclbcg.hardsepBCG.fits \
    ${path_wavelets}${dir}${run}${n:0:-5}.results.icl.hardsep.fits \
    ${path_wavelets}${dir}${run}${n:0:-5}.results.uicl.hardsep.fits \
    ${path_wavelets}${dir}${run}${n:0:-5}.results.gal.hardsep.fits \
    ${path_wavelets}${dir}${run}${n:0:-5}.results.residuals.fits \
    ${path_wavelets}${dir}${run}${n:0:-5}.results.rim.fits \
    -lock scalelimits yes \
    -lock scale yes \
    -lock frame image \
    -lock colorbar yes \
    -scale mode 99.5
  done
done
