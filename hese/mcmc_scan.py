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

import chainer_plot


BESTFIT = [1, 0, 0]
SIGMA = 0.2
NUFIT = False


def test_uni(x):
    """Test the unitarity of a matrix."""
    print 'Unitarity test:\n{0}'.format(abs(np.dot(x, x.conj().T)))


def angles_to_u(angles):
    s12_2, c13_4, s23_2, dcp = angles
    dcp = np.complex128(dcp)

    c13_2 = np.sqrt(c13_4)

    c12 = np.sqrt(1. - s12_2)
    s12 = np.sqrt(s12_2)
    c13 = np.sqrt(c13_2)
    s13 = np.sqrt(1. - c13_2)
    c23 = np.sqrt(1. - s23_2)
    s23 = np.sqrt(s23_2)

    phase = np.exp(-1j*dcp)
    p1 = np.array([[1   , 0   , 0]         , [0    , c23 , s23] , [0          , -s23 , c23]] , dtype=np.complex128)
    p2 = np.array([[c13 , 0   , s13*phase] , [0    , 1   , 0]   , [-s13*phase , 0    , c13]] , dtype=np.complex128)
    p3 = np.array([[c12 , s12 , 0]         , [-s12 , c12 , 0]   , [0          , 0    , 1]]   , dtype=np.complex128)

    u = np.dot(np.dot(p1, p2), p3)
    return u


def u_to_fr(initial_fr, matrix):
    """Compute the observed flavour ratio assuming decoherence."""
    composition = np.einsum(
        'ai, bi, a -> b', abs(matrix)**2, abs(matrix)**2, initial_fr
    )
    ratio = composition / np.sum(initial_fr)
    return ratio


def triangle_llh(theta):
    """-Log likelihood function for a given theta."""
    fr1, fr2 = theta[-2:]
    fr3 = 1.0 - (fr1 + fr2)

    u = angles_to_u(theta[:-2])
    fr = u_to_fr((fr1, fr2, fr3), u)
    fr_bf = BESTFIT
    cov_fr = np.identity(3) * SIGMA
    return -np.log10(multivariate_normal.pdf(fr, mean=fr_bf, cov=cov_fr))


def lnprior(theta):
    """Priors on theta."""
    s12_2, c13_4, s23_2, dcp, fr1, fr2 = theta

    fr3 = 1.0 - (fr1 + fr2)

    # Flavour ratio bounds
    if 0. <= fr1 <= 1.0 and 0. <= fr2 <= 1.0 and 0. <= fr3 <= 1.0:
        pass
    else: return -np.inf

    # Mixing angle bounds
    if 0. <= s12_2 <= 1. and 0. <= c13_4 <= 1. and 0. <= s23_2 <= 1. \
       and 0. <= dcp <= 2*np.pi:
        pass
    else: return -np.inf

    if NUFIT:
        u = angles_to_u(theta[:-2])
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
    print 'BESTFIT = {0}'.format(BESTFIT)
    print 'SIGMA = {0}'.format(SIGMA)
    print 'NUFIT = {0}'.format(NUFIT)

    ndim = 6
    nwalkers = args.nwalkers
    ntemps = 1
    burnin = args.burnin
    betas = np.array([1e0, 1e-1, 1e-2, 1e-3, 1e-4])
    if not NUFIT:
        p0_base = [0.5, 0.5, 0.5, np.pi, 0.5, 0.5]
        p0_std = [0.2] * ndim
    else:
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
    print sampler.acceptance_fraction
    print np.sum(sampler.acceptance_fraction)
    print np.unique(samples[:,0]).shape

    outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:03d}'.format(
        int(BESTFIT[0]*100), int(BESTFIT[1]*100), int(BESTFIT[2]*100), int(SIGMA*100)
    )
    if NUFIT:
        outfile += '_nufit'
    np.save(outfile+'.npy', samples)

    print "Making triangle plots"
    chainer_plot.plot(
        infile=outfile+'.npy',
        angles=False,
        nufit=NUFIT,
        outfile=outfile+'.pdf',
        bestfit_ratio=BESTFIT,
        sigma_ratio=SIGMA
    )
    chainer_plot.plot(
        infile=outfile+'.npy',
        angles=True,
        nufit=NUFIT,
        outfile=outfile+'_angles.pdf',
        bestfit_ratio=BESTFIT,
        sigma_ratio=SIGMA
    )
    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
