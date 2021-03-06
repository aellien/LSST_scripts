import dawis as d
import os
import sys
import glob
import numpy as np

if __name__ == '__main__':

    # Paths, lists & variables
    path_wavelets = '/n03data/ellien/LSST_ICL/wavelets/out2/'
    dirl = ['HorizonAGN', 'Hydrangea', 'Magneticum', 'TNG-100']

    lvl_sep_big = 5
    monomodality = False
    bspl = 1 / 16. * np.array([ 1, 4, 6, 4, 1 ])

    for dir in dirl:

        # Read Atoms
        rimpath = os.path.join(path_wavelets, dir, 'run1', '*.ol.*')
        grimpath = glob.glob(rimpath)
        grimpath.sort()

        for i, oln in enumerate( grimpath ):

            # Read iteration
            try:
                ol = d.read_objects_from_pickle( oln )
                itl = d.read_interscale_trees_from_pickle(  oln[:-12] + 'itl.it%03d.pkl' %(i+1) )
                wdc = d.wavelet_datacube.from_fits( oln[:-12] + 'wdc.it%03d.fits' %(i+1) )
                ldc = d.label_datacube.from_fits( oln[:-12] + 'ldc.it%03d.fits' %(i+1) )
            except:
                continue

            # If monomodality need to be done again
            for it in itl:
                if it.interscale_maximum.level >= lvl_sep_big :
                    if monomodality == True:
                        it, ldc = enforce_monomodality( it, wdc, ldc )

            # Read objects and interscale trees
            for j, ( o, it ) in enumerate( zip( ol, itl ) ):

                wt, st, lt = it.tubes( wdc, ldc )

                wr = wt - d.atrous( o.image, n_levels = wt.shape[2], \
                                            filter = bspl, \
                                            conditions = 'prolongation' )[1] * st

                norm_wr = np.sqrt( np.sum( wr**2 ) )
                sum_wr = np.sum( wr )

                o.norm_wr = norm_wr
                o.sum_wr = sum_wr

                print( np.sum(o.image), o.norm_wr, o.sum_wr )

            d.write_objects_to_pickle(ol, oln, overwrite = True)
