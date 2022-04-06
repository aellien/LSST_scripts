import dawis as d
import numpy as np
from astropy.io import fits
import os

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/ellien/LSST_ICL/simulations/out2/HorizonAGN'
    path_scripts = '/home/ellien/LSST_ICL/scripts'
    path_plots = '/home/ellien/LSST_ICL/plots/out2'
    path_wavelets = '/home/ellien/LSST_ICL/wavelets/out2/HorizonAGN/run2'

    n_levels = 11
    n_softhard_icl = 5
    n_hard_icl = 6

    bspl = 1 / 16. * np.array([ 1, 4, 6, 4, 1 ])

    nfl = [ '00761_0000001_0.05_xy_r_4Mpc_mu30.rebin.norm.results.icl.hardsep.fits', \
            '00761_0000001_0.05_yz_r_4Mpc_mu30.rebin.norm.results.icl.hardsep.fits' ]

    for nf in nfl:

        hdu = fits.open( os.path.join( path_wavelets, nf ))
        icl = hdu[0].data
        header = hdu[0].header

        cdc, wdc = d.bspl_atrous(icl, n_levels, header, 'prolongations')
        wdc.array[ :, :, :n_softhard_icl ] = 0.

        s = np.zeros( wdc.array.shape )
        s[ np.where( wdc.array > 0. )] = 1.
        sdc = d.datacube( s )

        ldc = d.label_regions( sdc )

        rl = d.make_regions_full_props( wdc, ldc )

        it = d.interscale_tree( rl[-1], rl, wdc, ldc )

        rec = it.CG_minimization( wdc, ldc, bspl, \
                            synthesis_operator = 'ADJOINT', step_size = 'FR', \
                            sigma_flux = 1E-3, max_iter = 200, verbose = True )

        hduo = fits.PrimaryHDU(rec, header)
        hduo.writeto( os.path.join( path_wavelets, ''.join((nf[:-4], 'post.fits')) ) )
