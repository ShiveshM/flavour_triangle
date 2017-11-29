#! /usr/bin/env python
"""
From a gaussian likelihood run an MCMC scan to find the posteriors
"""

from __future__ import absolute_import, division

import sys
sys.path.append('/users/mandalia/Documents/flavour_triangle/hese/mm_closure')

import argparse
import multiprocessing

import numpy as np
from scipy.stats import multivariate_normal

import emcee
import tqdm

import chainer_plot


FLAT = False
INJECTED = [2, 1]
BESTFIT = [1, 1]
SIGMA = 0.01
HYPO = [2, 1]
NUFIT = False


def test_uni(x):
    """Test the unitarity of a matrix."""
    print 'Unitarity test:\n{0}'.format(abs(np.dot(x, x.conj().T)))


def angles_to_u(angles):
    theta, delta = angles

    ct = np.cos(theta)
    st = np.sin(theta)

    u = np.array([[ct, -st*np.exp(-1j*delta)], [st*np.exp(1j*delta), ct]], dtype=np.complex128)
    return u


def u_to_fr(initial_fr, matrix):
    """Compute the observed flavour ratio assuming decoherence."""
    # TODO(shivesh): energy dependence
    composition = np.einsum(
        'ai, bi, a -> b', abs(matrix)**2, abs(matrix)**2, initial_fr
    )
    ratio = composition / np.sum(initial_fr)
    return ratio


def triangle_llh(angles):
    """-Log likelihood function for a given angles."""
    fr1, fr2 = HYPO

    u = angles_to_u(angles)
    fr = u_to_fr((fr1, fr2), u)
    fr_bf = BESTFIT
    cov_fr = np.identity(2) * SIGMA
    if FLAT:
        return 10.
    else:
        return np.log(multivariate_normal.pdf(fr, mean=fr_bf, cov=cov_fr))


def lnprior(angles):
    """Priors on angles."""
    theta, delta = angles
    # Mixing angle bounds
    if 0. <= theta <= 2*np.pi and 0. <= delta <= 2*np.pi:
        pass
    else: return -np.inf

    return 0.


def lnprob(angles):
    """Prob function for mcmc."""
    lp = lnprior(angles)
    if not np.isfinite(lp):
        return -np.inf
    return lp + triangle_llh(angles)


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--sigma-ratio', type=float, default=0.01,
        help='Set the 1 sigma for the flavour ratio'
    )
    parser.add_argument(
        '--injected-ratio', type=int, nargs=2, default=[2, 1],
        help='Set the injected initial (source) flavour ratio'
    )
    parser.add_argument(
        '--hypothesis-ratio', type=int, nargs=2, default=[2, 1],
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
        '--nsteps', type=int, default=2000,
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
        '--outfile', type=str, default='./untitled',
        help='Path to output chains'
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    np.random.seed(args.seed)

    global INJECTED
    global BESTFIT
    global HYPO
    global SIGMA
    global NUFIT
    INJECTED = np.array(args.injected_ratio) / float(np.sum(args.injected_ratio))
    HYPO = np.array(args.hypothesis_ratio) / float(np.sum(args.hypothesis_ratio))
    SIGMA = args.sigma_ratio
    if args.nufit.lower() == 'true':
        NUFIT = True
    elif args.nufit.lower() == 'false':
        NUFIT = False
    else:
        raise ValueError
    print 'HYPO = {0}'.format(HYPO)
    print 'INJECTED = {0}'.format(INJECTED)
    print 'SIGMA = {0}'.format(SIGMA)
    print 'NUFIT = {0}'.format(NUFIT)

    ANGLES = (1/4.*np.pi, 1.*np.pi)
    nufit_u = angles_to_u((ANGLES))
    print abs(nufit_u)
    BESTFIT = u_to_fr(INJECTED, nufit_u)

    ndim = 2
    nwalkers = args.nwalkers
    ntemps = 1
    burnin = args.burnin
    betas = np.array([1e0, 1e-1, 1e-2, 1e-3, 1e-4])
    if not NUFIT:
        p0_base = [np.pi, np.pi]
        p0_std = [0.2] * ndim
    else:
        assert 0
        p0_base = [0.306, 0.958, 0.441, 2*np.pi, 0.5, 0.5]
        p0_std = [0.005, 0.001, 0.01, 0.5, 0.2, 0.2]

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
    print 'acceptance fraction', sampler.acceptance_fraction
    print 'sum of acceptance fraction', np.sum(sampler.acceptance_fraction)
    print 'np.unique(samples[:,0]).shape', np.unique(samples[:,0]).shape

    print 'autocorrelation', sampler.acor

    outfile = args.outfile+'_{0:2f}_{1:2f}_{2:03d}_{3:03d}_{4:04d}'.format(
        abs(nufit_u[0][0]), abs(nufit_u[0][1]), int(INJECTED[0]*100), int(INJECTED[1]*100), int(SIGMA*1000)
    )
    if NUFIT:
        outfile += '_nufit'
    if FLAT:
        outfile += '_flat'
    # np.save(outfile+'.npy', samples)

    print "Making triangle plots"
    chainer_plot.plot(
        infile=outfile+'.npy',
        angles=True,
        outfile=outfile+'_angles.pdf',
        injected_ratio=INJECTED,
        hypothesis_ratio=HYPO,
        sigma_ratio=SIGMA
    )
    chainer_plot.plot(
        infile=outfile+'.npy',
        angles=False,
        outfile=outfile+'.pdf',
        injected_ratio=INJECTED,
        hypothesis_ratio=HYPO,
        sigma_ratio=SIGMA
    )
    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
