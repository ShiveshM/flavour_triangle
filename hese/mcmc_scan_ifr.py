#! /usr/bin/env python
"""
From a gaussian likelihood run an MCMC scan to find the posteriors
"""

from __future__ import absolute_import, division

import sys
sys.path.append('/users/mandalia/Documents/flavour_triangle/hese/')

import argparse
import multiprocessing

import numpy as np
from scipy.stats import multivariate_normal

import emcee
import tqdm

from mcmc_scan import test_uni, angles_to_u, u_to_fr
import chainer_ifr


BESTFIT = [1, 0, 0]
SIGMA = 0.2
HYPO = [1, 0, 0]
NUFIT = False


def triangle_llh(theta):
    """-Log likelihood function for a given theta."""
    fr1, fr2, fr3 = HYPO

    u = angles_to_u(theta)
    fr = u_to_fr((fr1, fr2, fr3), u)
    fr_bf = BESTFIT
    cov_fr = np.identity(3) * SIGMA
    llh = np.log10(multivariate_normal.pdf(fr, mean=fr_bf, cov=cov_fr))
    return -llh


def lnprior(theta):
    """Priors on theta."""
    s12_2, c13_4, s23_2, dcp = theta
    fr1, fr2, fr3 = HYPO

    # Mixing angle bounds
    if 0. <= s12_2 <= 1. and 0. <= c13_4 <= 1. and 0. <= s23_2 <= 1. \
       and 0. <= dcp <= 2*np.pi:
        pass
    else: return -np.inf

    if NUFIT:
        u = angles_to_u(theta)
        a_u = abs(u)
        ue1 = a_u[0][0]
        ue2 = a_u[0][1]
        ue3 = a_u[0][2]
        um1 = a_u[1][0]
        um2 = a_u[1][1]
        um3 = a_u[1][2]
        ut1 = a_u[2][0]
        ut2 = a_u[2][1]
        ut3 = a_u[2][2]

        # test_uni(u)
        # print a_u

        # Mixing elements 3sigma bound from nufit
        if 0.800 <= ue1 <= 0.844 and 0.515 <= ue2 <= 0.581 and 0.139 <= ue3 <= 0.155 \
        and 0.229 <= um1 <= 0.516 and 0.438 <= um2 <= 0.699 and 0.614 <= um3 <= 0.790 \
        and 0.249 <= ut1 <= 0.528 and 0.462 <= ut2 <= 0.715 and 0.595 <= ut3 <= 0.776:
            pass
        else:
            return -np.inf

    return 0.


def lnprob(theta):
    """Prob function for mcmc."""
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + triangle_llh(theta)


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--bestfit-ratio', type=int, nargs=3, default=[1, 0, 0],
        help='Set the bestfit flavour ratio'
    )
    parser.add_argument(
        '--sigma-ratio', type=float, default=0.2,
        help='Set the 1 sigma for the flavour ratio'
    )
    parser.add_argument(
        '--hypothesis-ratio', type=int, nargs=3, default=[1, 0, 0],
        help='Set the hypothesis for the flavour ratio'
    )
    parser.add_argument(
        '--burnin', type=int, default=100,
        help='Amount to burnin'
    )
    parser.add_argument(
        '--nwalkers', type=int, default=100,
        help='Number of walkers'
    )
    parser.add_argument(
        '--nsteps', type=int, default=1000,
        help='Number of steps to run'
    )
    parser.add_argument(
        '--nufit', type=str, default='False',
        help='Include NuFit priors'
    )
    parser.add_argument(
        '--seed', type=int, default=99,
        help='Set the random seed value'
    )
    parser.add_argument(
        '--outfile', type=str, default='./untitled.npy',
        help='Path to output chains'
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    np.random.seed(args.seed)

    global BESTFIT
    global SIGMA
    global NUFIT
    BESTFIT = np.array(args.bestfit_ratio) / float(np.sum(args.bestfit_ratio))
    SIGMA = args.sigma_ratio
    if args.nufit.lower() == 'true':
        NUFIT = True
    elif args.nufit.lower() == 'false':
        NUFIT = False
    else:
        raise ValueError

    HYPO = np.array(args.hypothesis_ratio) / float(np.sum(args.hypothesis_ratio))

    print 'HYPO = {0}'.format(HYPO)
    print 'BESTFIT = {0}'.format(BESTFIT)
    print 'SIGMA = {0}'.format(SIGMA)
    print 'NUFIT = {0}'.format(NUFIT)

    ndim = 4
    nwalkers = args.nwalkers
    ntemps = 1
    burnin = args.burnin
    betas = np.array([1e0, 1e-1, 1e-2, 1e-3, 1e-4])
    if not NUFIT:
        p0_base = [0.5, 0.5, 0.5, np.pi]
        p0_std = [0.2] * ndim
    else:
        p0_base = [0.306, 0.958, 0.441, 2*np.pi]
        p0_std = [0.005, 0.001, 0.01, 0.5]

    p0 = np.random.normal(p0_base, p0_std, size=[ntemps, nwalkers, ndim])
    print map(lnprior, p0[0])

    # threads = multiprocessing.cpu_count()
    threads = 1
    sampler = emcee.PTSampler(
        ntemps, nwalkers, ndim, triangle_llh, lnprior, threads=threads
    )

    print "Running burn-in"
    for result in tqdm.tqdm(sampler.sample(p0, iterations=burnin), total=burnin):
        pos, prob, state = result
    sampler.reset()
    print "Finished burn-in"

    nsteps = args.nsteps

    print "Running"
    for _ in tqdm.tqdm(sampler.sample(pos, iterations=nsteps), total=nsteps):
        pass
    print "Finished"

    samples = sampler.chain[0, :, :, :].reshape((-1, ndim))
    print sampler.acceptance_fraction
    print np.sum(sampler.acceptance_fraction)
    print np.unique(samples[:,0]).shape

    outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:03d}_ifr_{4:03d}_{5:03d}_{6:03d}'.format(
        int(BESTFIT[0]*100), int(BESTFIT[1]*100), int(BESTFIT[2]*100), int(SIGMA*100),
        int(HYPO[0]*100), int(HYPO[1]*100), int(HYPO[2]*100)
    )
    if NUFIT:
        outfile += '_nufit'
    np.save(outfile+'.npy', samples)

    print "Making triangle plots"
    chainer_ifr.plot(
        infile=outfile+'.npy',
        angles=False,
        nufit=NUFIT,
        outfile=outfile+'.pdf',
        hypothesis_ratio=HYPO,
        bestfit_ratio=BESTFIT,
        sigma_ratio=SIGMA
    )
    chainer_ifr.plot(
        infile=outfile+'.npy',
        angles=True,
        nufit=NUFIT,
        outfile=outfile+'_angles.pdf',
        hypothesis_ratio=HYPO,
        bestfit_ratio=BESTFIT,
        sigma_ratio=SIGMA
    )
    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
