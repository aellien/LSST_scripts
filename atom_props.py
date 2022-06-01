# Last modif: 08/2021
# Author: Amael Ellien
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Modules
import dawis as d
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import pandas as pd
from astropy.io import fits
from astropy.visualization import LinearStretch, LogStretch
from astropy.visualization import ZScaleInterval, MinMaxInterval
from astropy.visualization import ImageNormalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import SymLogNorm
from scipy.stats import sigmaclip
from skimage.measure import label, regionprops
from sklearn.utils import resample

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def average_size_atom( ol, n_levels ):

    sizes = np.zeros( (n_levels, 4) )

    for o in ol:
        x_min, y_min, x_max, y_max = o.bbox
        xs = x_max - x_min
        ys = y_max - y_min
        sizes[o.level, 0] += xs
        sizes[o.level, 1] += 1
        sizes[o.level, 2] += ys
        sizes[o.level, 3] += 1

    print(nf)
    data = []
    for i in range(n_levels):
        data.append( [ 2**i, sizes[i, 0] / sizes[i, 1], sizes[i, 2] / sizes[i, 3]] )

    df = pd.DataFrame( data, columns = [ 'z', '<sx>', '<sy>'] )

    return df


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/n03data/ellien/LSST_ICL/simulations/out2/'
    path_scripts = '/home/ellien/LSST_ICL/LSST_scripts'
    path_wavelets = '/n03data/ellien/LSST_ICL/wavelets/out2/'
    path_plots = '/n03data/ellien/LSST_ICL/plots'

    dirl = ['HorizonAGN', 'Hydrangea', 'Magneticum', 'TNG-100']

    gamma = 0.5
    n_levels = 11
    lvl_sep_big = 5
    n_hard_icl = 5
    pixscale = 1.6 # ''/pixel
    physscale = 1 # kpc/''
    rc = 20 # pixels, distance to center to be classified as gal
    ricl = 650 # pixels, distance to center to be classified as ICL

    plt.figure()

    for dir in dirl:

        image_dir = os.path.join( path_data, dir )
        image_files = glob.glob(image_dir+'/*norm.fits')

        for nf in image_files:

            print('\n%s'%nf)

            # Read atoms
            nf = nf[:-4]
            split = nf.split('/')
            nf = split[-1]
            nfp = os.path.join( path_wavelets, dir, 'run1', nf )
            opath = nfp + '*ol*pkl'
            opathl = glob.glob(opath)
            opathl.sort()
            itpath = nfp + '*itl*pkl'
            itpathl = glob.glob(itpath)
            itpathl.sort()
            tol = []
            titl = []

            for i, ( op, itlp ) in enumerate( zip( opathl, itpathl )):

                print('read iteration %d' %(i), end ='\r')

                ol = d.read_objects_from_pickle( op )
                itl = d.read_interscale_trees_from_pickle( itlp )

                for j, o in enumerate(ol):

                    tol.append(o)
                    titl.append(itl[j])

            print(len(tol))

            df = average_size_atom( tol, n_levels )
            plt.plot( df['z'], df['sx'], color = 'blue', alpha = 0.5 )
            plt.plot( df['z'], df['sx'], color = 'red', alpha = 0.5 )

    plt.savefig(os.path.join(path_plots, 'average_size_vs_z.pdf'), format = 'pdf')
