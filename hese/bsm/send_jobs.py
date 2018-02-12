#! /usr/bin/env python

import numpy as np
import os

a_fr = (1, 2, 0)
b_fr = (1, 0, 0)
c_fr = (0, 1, 0)
d_fr = (0, 0, 1)
e_fr = (1, 1, 1)
f_fr = (2, 1, 0)
g_fr = (1, 1, 0)

full_scan_mfr = [
    (1, 1, 1), (1, 1, 0)
]

fix_sfr_mfr = [
    (1, 1, 1, 1, 0, 0),
    (1, 1, 1, 0, 1, 0),
    (1, 1, 1, 0, 0, 1),
    (1, 1, 1, 1, 2, 0),
    (1, 1, 0, 0, 1, 0),
    (1, 1, 0, 1, 2, 0),
    (1, 1, 0, 1, 0, 0),
    (1, 0, 0, 1, 0, 0),
    (0, 1, 0, 0, 1, 0),
    (1, 2, 0, 0, 1, 0),
    (1, 2, 0, 1, 2, 0)
]

sigmas = ['0.1', '0.01']
dimensions = [3, 6]
spectral_index = ['-2']
binning = [1e4, 1e7, 20]
flat = False
burnin = 500
nwalkers = 200
nsteps = 5000
scales = "1E-20 1E-30"

script = '/users/mandalia/Documents/flavour_triangle/hese/bsm/wrap.sh'

job_number = 1
for dim in dimensions:
    print 'dimension', dim
    for s_i in spectral_index:
        print 'spectral index {0}'.format(s_i)

        outchain_head = '/data/icecube/mandalia/flavour_ratio/data/DIM{0}/SI_{1}'.format(dim, s_i)

        for sig in sigmas:
            print 'sigma', sig
            for frs in fix_sfr_mfr:
                print frs
                outchains = outchain_head + '/fix_ifr/{0}/mcmc_chain'.format(str(sig).replace('.', '_'))
                command = ''
                command += 'qsub -cwd -V -q SL6 '
                command += script
                command += ' {0}'.format(frs[0])
                command += ' {0}'.format(frs[1])
                command += ' {0}'.format(frs[2])
                command += ' {0}'.format(sig)
                command += ' {0}'.format('True')
                command += ' {0}'.format(frs[3])
                command += ' {0}'.format(frs[4])
                command += ' {0}'.format(frs[5])
                command += ' {0}'.format('False')
                command += ' {0}'.format(0)
                command += ' {0}'.format(dim)
                command += ' {0}'.format(s_i)
                command += ' {0}'.format(flat)
                command += ' {0}'.format(burnin)
                command += ' {0}'.format(nwalkers)
                command += ' {0}'.format(nsteps)
                command += ' {0}'.format(outchains)
                command += ' {0}'.format('False')
                command += ' {0}'.format(binning[0])
                command += ' {0}'.format(binning[1])
                command += ' {0}'.format(binning[2])
                os.system(command)

            for frs in full_scan_mfr:
                print frs
                outchains = outchain_head + '/full_scan/{0}/mcmc_chain'.format(str(sig).replace('.', '_'))
                command = ''
                command += 'qsub -cwd -V -q SL6 '
                command += script
                command += ' {0}'.format(frs[0])
                command += ' {0}'.format(frs[1])
                command += ' {0}'.format(frs[2])
                command += ' {0}'.format(sig)
                command += ' {0}'.format('False')
                command += ' {0}'.format(0)
                command += ' {0}'.format(0)
                command += ' {0}'.format(0)
                command += ' {0}'.format('False')
                command += ' {0}'.format(0)
                command += ' {0}'.format(dim)
                command += ' {0}'.format(s_i)
                command += ' {0}'.format(flat)
                command += ' {0}'.format(burnin)
                command += ' {0}'.format(nwalkers)
                command += ' {0}'.format(nsteps)
                command += ' {0}'.format(outchains)
                command += ' {0}'.format('False')
                command += ' {0}'.format(binning[0])
                command += ' {0}'.format(binning[1])
                command += ' {0}'.format(binning[2])
                os.system(command)
