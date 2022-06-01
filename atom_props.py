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
from scipy.stats import sigmaclip, median_absolute_deviation
from skimage.measure import label, regionprops
from sklearn.utils import resample

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def average_size_atom( ol, n_levels ):

    sx = list( np.zeros( ( n_levels, 1 ) ))
    sy = list( np.zeros( ( n_levels, 1 ) ))

    for o in ol:
        x_min, y_min, x_max, y_max = o.bbox
        xs = x_max - x_min
        ys = y_max - y_min
        sx[o.level] = np.append( sx[o.level], xs )
        sy[o.level] = np.append( sy[o.level], ys )

    data = []
    for i in range(n_levels):
        data.append( [ i, 2**i, np.mean(sx[i]), np.mean(sy[i]), np.std(sx[i]), np.std(sy[i]), np.median(sx[i]), np.median(sy[i]), median_absolute_deviation(sx[i]), median_absolute_deviation([sy[i]]) ] )

    df = pd.DataFrame( data, columns = [ 'z', 'sz', '<sx>', '<sy>', 'std(sx)', 'std(sy)', 'med(sx)', 'med(sy)', 'mad(sx)', 'mad(sy)'] )

    return df


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/n03data/ellien/LSST_ICL/simulations/out2/'
    path_scripts = '/home/ellien/LSST_ICL/LSST_scripts'
    path_wavelets = '/n03data/ellien/LSST_ICL/wavelets/out2/'
    path_plots = '/home/ellien/LSST_ICL/plots'

    dirl = ['HorizonAGN', 'Hydrangea', 'Magneticum', 'TNG-100']

    gamma = 0.5
    n_levels = 11
    lvl_sep_big = 5
    n_hard_icl = 5
    pixscale = 1.6 # ''/pixel
    physscale = 1 # kpc/''
    rc = 20 # pixels, distance to center to be classified as gal
    ricl = 650 # pixels, distance to center to be classified as ICL

    fig, ax_tot = plt.subplots( 1 )

    dfl = []

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

            figo, ax = plt.subplots( 1 )

            for i, ( op, itlp ) in enumerate( zip( opathl, itpathl )):

                print('read iteration %d' %(i), end ='\r')

                ol = d.read_objects_from_pickle( op )
                itl = d.read_interscale_trees_from_pickle( itlp )

                for j, o in enumerate(ol):

                    tol.append(o)
                    titl.append(itl[j])

            print('\n%d sources in total.' %len(tol))

            df = average_size_atom( tol, n_levels )
            dfl.append(df)

            ax.plot( df['z'], df['<sx>'], color = 'blue', alpha = 0.5 )
            ax.plot( df['z'], df['<sx>'], color = 'blue', alpha = 0.5 )
            ax.plot( df['z'], df['med(sx)'], color = 'red', alpha = 0.5 )
            ax.plot( df['z'], df['med(sy)'], color = 'red', alpha = 0.5 )
            figo.savefig( os.path.join( path_wavelets, dir, 'run1', nf + 'average_size_vs_z.pdf' ), format = 'pdf')
            plt.close()

    sxl = list( np.zeros( ( n_levels, 1 ) ))
    syl = list( np.zeros( ( n_levels, 1 ) ))
    mxl = list( np.zeros( ( n_levels, 1 ) ))
    myl = list( np.zeros( ( n_levels, 1 ) ))


    for i in range(n_levels):
        for df in dfl:
            sxl[i] = np.append( sxl[i], df['<sx>'][i] )
            syl[i] = np.append( syl[i], df['<sy>'][i] )
            mxl[i] = np.append( sxl[i], df['med(sx)'][i] )
            myl[i] = np.append( syl[i], df['med(sy)'][i] )

    ax_tot.errorbar( df['z'], [ np.mean(x) for x in sxl ], yerr = [ np.std(x) for x in sxl ], color = 'blue', alpha = 0.5 )
    ax_tot.errorbar( df['z'], [ np.mean(y) for y in syl ], yerr = [ np.std(y) for y in syl ], color = 'red', alpha = 0.5 )
    ax_tot.errorbar( df['z'], [ np.mean(x) for x in mxl ], yerr = [ np.std(x) for x in mxl ], color = 'green', alpha = 0.5 )
    ax_tot.errorbar( df['z'], [ np.mean(y) for y in myl ], yerr = [ np.std(y) for y in myl ], color = 'pink', alpha = 0.5 )
    fig.savefig(os.path.join(path_plots, 'average_size_vs_z.pdf'), format = 'pdf')

    plt.close()
