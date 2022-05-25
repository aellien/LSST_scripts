import dawis as d
import os
import sys
import glob

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/ellien/LSST_ICL/simulations/out2/'
    path_scripts = '/home/ellien/LSST_ICL/scripts'
    path_wavelets = '/home/ellien/LSST_ICL/wavelets/out2/'

    dirl = ['Hydrangea']

    framerate = 3

    for dir in dirl:

        image_dir = os.path.join( path_data, dir )
        image_files = glob.glob(image_dir+'/*.fits')

        for nf in image_files:

            split = nf.split('/')
            nf = split[-1]
            nfp = os.path.join( path_wavelets, dir, 'run1', nf[:-5] )
            d.make_gif( framerate, nfp )
