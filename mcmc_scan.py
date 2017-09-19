#! /usr/bin/env python
"""
From a gaussian likelihood run an MCMC scan to find the posteriors
"""

from __future__ import absolute_import, division

import argparse
import multiprocessing

import numpy as np
from scipy.stats import multivariate_normal

import emcee
import tqdm


BESTFIT = [0.5, 0.5, 0.0]
SIGMA = 0.2


def triangle_llh(theta):
    """-Log likelihood function for a given theta."""
    fr = theta[-3:]
    fr_bf = BESTFIT
    SIGMA = 0.2
    cov_fr = np.identity(3) * SIGMA
    return -np.log10(multivariate_normal.pdf(fr, mean=fr_bf, cov=cov_fr))


def lnprior(theta):
    """Priors on theta."""
    ue1, ue2, ue3, um1, um2, um3, ut1, ut2, ut3, fr1, fr2, fr3 = theta
    
    allow = True
    # Flavour ratio bounds
    if 0. <= fr1 <= 1.0 and 0. <= fr2 <= 1.0 and 0. <= fr3 <= 1.0:
        pass
    else: allow = False

    # mixing elements 3sigma bound using nufit
    if 0.800 <= ue1 <= 0.844 and 0.515 <= ue2 <= 0.581 and 0.139 <= ue3 <= 0.155 \
    and 0.229 <= um1 <= 0.516 and 0.438 <= um2 <= 0.699 and 0.614 <= um3 <= 0.790 \
    and 0.249 <= ut1 <= 0.528 and 0.462 <= ut2 <= 0.715 and 0.595 <= ut3 <= 0.776:
        pass
    else: allow = False
        
    # TODO(shivesh): enforce sum of flavour ratios? enforce unitarity?
    if allow: return 0.
    else: return -np.inf


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
        '--sigma-ratio', type=float, default=0.2,
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
        '--nsteps', type=int, default=1000,
        help='Number of steps to run'
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

    BESTFIT = np.array(args.bestfit_ratio) / float(np.sum(args.bestfit_ratio))
    SIGMA = args.sigma_ratio

    ndim = 12
    nwalkers = args.nwalkers
    ntemps = 1
    burnin = args.burnin
    betas = np.array([1e0, 1e-1, 1e-2, 1e-3, 1e-4])
    p0_base = [0.82, 0.55, 0.14, 0.40, 0.50, 0.65, 0.40, 0.60, 0.65, 0.5, 0.5, 0.5]
    p0_std = [0.01, 0.01, 0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.1, 0.1, 0.1]

    p0 = np.random.normal(p0_base, p0_std, size=[ntemps, nwalkers, ndim])
    print map(lnprior, p0[0])

    thread = multiprocessing.cpu_count()
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
    print np.unique(samples[:,0]).shape

    np.save(args.outfile, samples)


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
