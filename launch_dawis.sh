#!/bin/bash
#
# 00761_0000001_0.05_xy_r_4Mpc_mu30.rebin.norm.fits
# 00761_0000001_0.05_yz_r_4Mpc_mu30.rebin.norm.fits
#

for file in 00761_0000001_0.05_xy_r_4Mpc_mu30.rebin.norm.fits 00761_0000001_0.05_yz_r_4Mpc_mu30.rebin.norm.fits
do
      echo "Launch Dawis on file ${file}"
      qsub qsub_dawis.sh -v ncl=${file:12:7},nf=${file}
      sleep 2
done
