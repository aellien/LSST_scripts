import os
import dawis
import sys
import shutil

indir = os.path.join( '/n03data/ellien/LSST_ICL/simulations/out4/', sys.argv[2] )
infile = sys.argv[1]
outdir = os.path.join( '/n03data/ellien/LSST_ICL/wavelets/out4/', sys.argv[2], 'run1' )
n_cpus = 4 # Number of CPUs
tau = 0.1   # Relative Threshold /!\ different des autres simus
gamma = 0.8   # Attenuation (CLEAN) factor
ceps = 1E-4    # Convergence value for epsilon
n_levels = 11    # Number of wavelet scales
min_span = 1    # Minimum of wavelet scales spanned by an interscale tree (must be >= 1)
max_span = 2    # Maximum number of wavelet scales spanned by an interscale tree
lvl_sep_big = 6     # Scale at wich mix_span & max_span are set to 1, and gamma to 1
extent_sep = 1E-10    # Ratio n_pix/vignet under which the Haar wavelet is used for restoration
lvl_sep_lin = -1     # Wavelet scale under which the Haar wavelet can be used for restoration
max_iter = 1500      # Maximum number of iterations
data_dump = True    # Write data at each iteration /!\ demands lot of space on hardware /!\
gif = True      # Make gifs of the run (need data_dump = True)
starting_level = 2 # Starting wavelet scale (this is the third scale - Python convention 0 1 2)
conditions = 'prolongation'
monomodality = False
resume = True
rm_gamma_for_big = True

shutil.copyfile( os.path.abspath(__file__), os.path.join( outdir, infile[:-4] + 'input.dawis.py' ) )

dawis.synthesis_by_analysis( indir = indir, infile = infile, outdir = outdir, n_cpus = n_cpus, n_levels = n_levels, \
                                    tau = tau, gamma = gamma, ceps = ceps, conditions = conditions, min_span = min_span, \
                                    max_span = max_span, lvl_sep_big = lvl_sep_big, monomodality = monomodality, \
                                    max_iter = max_iter, extent_sep = extent_sep, rm_gamma_for_big = rm_gamma_for_big, \
                                    data_dump = data_dump, gif = gif, resume = resume )
