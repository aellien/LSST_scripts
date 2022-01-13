import os
import dawis
import sys

indir = '/home/ellien/LSST_ICL/simulations/out1'
infile = sys.argv[1]
outdir = '/n03data/ellien/LSST_ICL/wavelets/out4'
n_cpus = 4 # Number of CPUs
tau = 0.1   # Relative Threshold  0.8
gamma = 0.2   # Attenuation (CLEAN) factor  0.2
ceps = 1E-5    # Convergence value for epsilon  1E-5
n_levels = 11    # Number of wavelet scales
min_span = 2    # Minimum of wavelet scales spanned by an interscale tree (must be >= 1)
max_span = 3    # Maximum number of wavelet scales spanned by an interscale tree
lvl_sep_big = 7     # Scale at wich mix_span & max_span are set to 1, and gamma to 1
extent_sep = 1.0    # Ratio n_pix/vignet under which the Haar wavelet is used for restoration
lvl_sep_lin = 2     # Wavelet scale under which the Haar wavelet can be used for restoration
max_iter = 1500      # Maximum number of iterations
data_dump = True    # Write data at each iteration /!\ demands lot of space on hardware /!\
gif = True      # Make gifs of the run (need data_dump = True)
starting_level = 2 # Starting wavelet scale (this is the third scale - Python convention 0 1 2)
conditions = 'prolongation'

dawis.synthesis_by_analysis( indir = indir, infile = infile, outdir = outdir, n_cpus = n_cpus, n_levels = n_levels, \
                                    tau = tau, gamma = gamma, ceps = ceps, conditions = conditions, min_span = min_span, \
                                    max_span = max_span, lvl_sep_big = lvl_sep_big, max_iter = max_iter, \
                                    data_dump = data_dump, gif = gif )
