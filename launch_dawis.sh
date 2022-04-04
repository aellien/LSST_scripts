#!/bin/bash
#
# HorizonAGN  Hydrangea  Magneticum  TNG-100

path="/n03data/ellien/LSST_ICL/simulations/out2/"

for dir in Hydrangea  Magneticum  TNG-100
do
    for file in $path/$dir/*norm.fits
    do
          echo "Launch Dawis on file ${file}"
          n=$(basename "$file")
          echo "qsub qsub_dawis.sh -v ncl=${n:0:21},nf=${n},dir=${dir}"
          sleep 2
    done
done
