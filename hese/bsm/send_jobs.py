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

sigmas = ['0.1', '0.01', '0.001', '0.0001', '0.00001']
dimensions = [3, 4, 5, 6]
energy = 100000
flat = False
burnin = 2000
nwalkers = 200
nsteps = 20000
scales = "1E-20 1E-30"

script = '/users/mandalia/Documents/flavour_triangle/hese/bsm/wrap.sh'

job_number = 1
for dim in dimensions:
    print 'dimension', dim
    outchain_head = '/data/mandalia/flavour_ratio/data/DIM{0}'.format(dim)

    for sig in sigmas:
        print 'sigma', sig
        for frs in fix_sfr_mfr:
            print frs
            outchains = outchain_head + '/fix_ifr/{0}/mcmc_chain'.format(str(sig).replace('.', '_'))
            command = ''
            # command += 'qsub -cwd -V '
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
            command += ' {0}'.format(energy)
            command += ' {0}'.format(flat)
            command += ' {0}'.format(burnin)
            command += ' {0}'.format(nwalkers)
            command += ' {0}'.format(nsteps)
            command += ' {0}'.format(outchains)
            command += ' {0}'.format('False')
            os.system(command)

        for frs in full_scan_mfr:
            print frs
            outchains = outchain_head + '/full_scan/{0}/mcmc_chain'.format(str(sig).replace('.', '_'))
            command = ''
            # command += 'qsub -cwd -V '
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
            command += ' {0}'.format(energy)
            command += ' {0}'.format(flat)
            command += ' {0}'.format(burnin)
            command += ' {0}'.format(nwalkers)
            command += ' {0}'.format(nsteps)
            command += ' {0}'.format(outchains)
            command += ' {0}'.format('False')
            os.system(command)
