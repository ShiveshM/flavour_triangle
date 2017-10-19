#! /usr/bin/env python
"""
From a gaussian likelihood run an MCMC scan to find the posteriors
"""

from __future__ import absolute_import, division

import sys
sys.path.append('/users/mandalia/Documents/flavour_triangle/hese/test_haar')

import argparse
import multiprocessing

import numpy as np
from scipy.stats import multivariate_normal

import emcee
import tqdm

import chainer_plot


BESTFIT = [1, 1, 1]
SIGMA = 0.01


def triangle_llh(theta):
    """-Log likelihood function for a given theta."""
    fr = theta

    fr_bf = BESTFIT
    cov_fr = np.identity(3) * SIGMA
    return np.log(multivariate_normal.pdf(fr, mean=fr_bf, cov=cov_fr))


def lnprior(theta):
    """Priors on theta."""
    fr1, fr2, fr3 = theta

    # Flavour ratio bounds
    if 0. <= fr1 <= 1.0 and 0. <= fr2 <= 1.0 and 0. <= fr3 <= 1.0:
        pass
    else: return -np.inf

    if fr1 + fr2 + fr3 <= 1.0:
        pass
    else: return -np.inf

    if abs(1. - (fr1 + fr2 + fr3)) < 0.01:
        pass
    else: return -np.inf

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
        '--bestfit-ratio', type=int, nargs=3, default=[1, 1, 1],
        help='Set the bestfit flavour ratio'
    )
    parser.add_argument(
        '--sigma-ratio', type=float, default=0.01,
        help='Set the 1 sigma for the flavour ratio'
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
        '--nsteps', type=int, default=2000,
        help='Number of steps to run'
    )
    parser.add_argument(
        '--seed', type=int, default=99,
        help='Set the random seed value'
    )
    parser.add_argument(
        '--outfile', type=str, default='./untitled',
        help='Path to output chains'
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    np.random.seed(args.seed)

    global BESTFIT
    global SIGMA
    BESTFIT = np.array(args.bestfit_ratio) / float(np.sum(args.bestfit_ratio))
    SIGMA = args.sigma_ratio
    print 'BESTFIT = {0}'.format(BESTFIT)
    print 'SIGMA = {0}'.format(SIGMA)

    ndim = 3
    nwalkers = args.nwalkers
    ntemps = 1
    burnin = args.burnin
    betas = np.array([1e0, 1e-1, 1e-2, 1e-3, 1e-4])
    p0_base = [0.33, 0.33, 0.33]
    p0_std = [0.1] * ndim

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

    outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:04d}'.format(
        int(BESTFIT[0]*100), int(BESTFIT[1]*100), int(BESTFIT[2]*100), int(SIGMA*1000)
    )
    np.save(outfile+'.npy', samples)

    print "Making triangle plots"
    chainer_plot.plot(
        infile=outfile+'.npy',
        angles=True,
        outfile=outfile+'_angles.pdf',
        bestfit_ratio=BESTFIT,
        sigma_ratio=SIGMA
    )
    # chainer_plot.plot(
    #     infile=outfile+'.npy',
    #     angles=False,
    #     outfile=outfile+'.pdf',
    #     bestfit_ratio=BESTFIT,
    #     sigma_ratio=SIGMA
    # )
    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
