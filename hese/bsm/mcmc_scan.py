#! /usr/bin/env python
"""
From a gaussian likelihood run an MCMC scan to find the posteriors
"""

from __future__ import absolute_import, division

import sys
sys.path.append('/users/mandalia/Documents/flavour_triangle/hese/bsm/')

import argparse
import multiprocessing

import numpy as np
from scipy.stats import multivariate_normal
from numpy import linalg as LA

import emcee
import tqdm

import chainer_plot


DIMENSION = 6
ENERGY = 100 # GeV
BESTFIT = [1, 1, 1]
SIGMA = 0.01
NUFIT = False
FLAT = True


def test_uni(x):
    """Test the unitarity of a matrix."""
    print 'Unitarity test:\n{0}'.format(abs(np.dot(x, x.conj().T)))

def angles_to_fr(angles):
    sphi4, c2psi = angles

    psi = (0.5)*np.arccos(c2psi)

    sphi2 = np.sqrt(sphi4)
    cphi2 = 1. - sphi2
    spsi2 = np.sin(psi)**2
    cspi2 = 1. - spsi2

    x = sphi2*cspi2
    y = sphi2*spsi2
    z = cphi2
    return x, y, z


def angles_to_u(angles):
    s12_2, c13_4, s23_2, dcp = angles
    dcp = np.complex128(dcp)

    c12_2 = 1. - s12_2
    c13_2 = np.sqrt(c13_4)
    s13_2 = 1. - c13_2
    c23_2 = 1. - s23_2

    # TODO(shivesh): negative sign?
    t12 = np.arcsin(np.sqrt(s12_2))
    t13 = np.arccos(np.sqrt(c13_2))
    t23 = np.arcsin(np.sqrt(s23_2))

    c12 = np.cos(t12)
    s12 = np.sin(t12)
    c13 = np.cos(t13)
    s13 = np.sin(t13)
    c23 = np.cos(t23)
    s23 = np.sin(t23)

    p1 = np.array([[1   , 0   , 0]                   , [0    , c23 , s23] , [0                   , -s23 , c23]] , dtype=np.complex128)
    p2 = np.array([[c13 , 0   , s13*np.exp(-1j*dcp)] , [0    , 1   , 0]   , [-s13*np.exp(1j*dcp) , 0    , c13]] , dtype=np.complex128)
    p3 = np.array([[c12 , s12 , 0]                   , [-s12 , c12 , 0]   , [0                   , 0    , 1]]   , dtype=np.complex128)

    u = np.dot(np.dot(p1, p2), p3)
    return u

NUFIT_U = angles_to_u((0.307, (1-0.2195)**2, 0.565, 3.97935))

def params_to_BSMu(theta):
    s12_2, c13_4, s23_2, dcp, sc1, sc2 = theta

    mass_matrix = np.array(
        [[0, 0, 0], [0, 7.40E-23, 0], [0, 0, 2.515E-21]]
    )
    sm_ham = (1./2*ENERGY)*np.dot(NUFIT_U, np.dot(mass_matrix**2, NUFIT_U.conj()))

    new_physics_u = angles_to_u((s12_2, c13_4, s23_2, dcp))
    scale_matrix = np.array(
        [[0, 0, 0], [0, sc1, 0], [0, 0, sc2]]
    )
    bsm_term = (ENERGY**(DIMENSION-3)) * np.dot(new_physics_u, np.dot(scale_matrix, new_physics_u.conj()))

    bsm_ham = sm_ham + bsm_term

    eg_values, eg_vector = LA.eig(bsm_ham)
    return eg_vector

def u_to_fr(initial_fr, matrix):
    """Compute the observed flavour ratio assuming decoherence."""
    # TODO(shivesh): energy dependence
    composition = np.einsum(
        'ai, bi, a -> b', abs(matrix)**2, abs(matrix)**2, initial_fr
    )
    ratio = composition / np.sum(initial_fr)
    return ratio


def triangle_llh(theta):
    """-Log likelihood function for a given theta."""
    fr1, fr2, fr3 = angles_to_fr(theta[-2:])

    u = params_to_BSMu(theta[:-2])
    fr = u_to_fr((fr1, fr2, fr3), u)
    fr_bf = BESTFIT
    cov_fr = np.identity(3) * SIGMA
    if FLAT:
        return 10.
    else:
        return np.log(multivariate_normal.pdf(fr, mean=fr_bf, cov=cov_fr))


def lnprior(theta):
    """Priors on theta."""
    s12_2, c13_4, s23_2, dcp, sc1, sc2, sphi4, c2psi = theta
    sc1, sc2 = np.power(10., sc1), np.power(10., sc2)

    # Flavour ratio bounds
    if 0. <= sphi4 <= 1.0 and -1.0 <= c2psi <= 1.0:
        pass
    else: return -np.inf

    # Mixing angle bounds
    if 0. <= s12_2 <= 1. and 0. <= c13_4 <= 1. and 0. <= s23_2 <= 1. \
       and 0. <= dcp <= 2*np.pi:
        pass
    else: return -np.inf

    # Scale bounds
    if 1e-40 <= sc1 <= 1e-30 and 1e-40 <= sc2 <= 1e-30:
        pass
    else: return -np.inf

    # if NUFIT:
    #     u = angles_to_u(theta[:-2])
    #     a_u = abs(u)
    #     ue1 = a_u[0][0]
    #     ue2 = a_u[0][1]
    #     ue3 = a_u[0][2]
    #     um1 = a_u[1][0]
    #     um2 = a_u[1][1]
    #     um3 = a_u[1][2]
    #     ut1 = a_u[2][0]
    #     ut2 = a_u[2][1]
    #     ut3 = a_u[2][2]

    #     # Mixing elements 3sigma bound from nufit
    #     if 0.800 <= ue1 <= 0.844 and 0.515 <= ue2 <= 0.581 and 0.139 <= ue3 <= 0.155 \
    #     and 0.229 <= um1 <= 0.516 and 0.438 <= um2 <= 0.699 and 0.614 <= um3 <= 0.790 \
    #     and 0.249 <= ut1 <= 0.528 and 0.462 <= ut2 <= 0.715 and 0.595 <= ut3 <= 0.776:
    #         pass
    #     else:
    #         return -np.inf

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
        assert 0
    elif args.nufit.lower() == 'false':
        NUFIT = False
    else:
        raise ValueError
    print 'BESTFIT = {0}'.format(BESTFIT)
    print 'SIGMA = {0}'.format(SIGMA)
    print 'NUFIT = {0}'.format(NUFIT)

    ndim = 8
    nwalkers = args.nwalkers
    ntemps = 1
    burnin = args.burnin
    betas = np.array([1e0, 1e-1, 1e-2, 1e-3, 1e-4])
    if not NUFIT:
        p0_base = [0.5, 0.5, 0.5, np.pi, -35, -35, 0.5, 0.0]
        p0_std = [0.2, 0.2, 0.2, 0.2, 5, 5, 0.2, 0.2]
    else:
        assert 0
        p0_base = [0.306, 0.958, 0.441, 2*np.pi, 0.5, 0.5]
        p0_std = [0.005, 0.001, 0.01, 0.5, 0.2, 0.2]

    print 'p0_base', p0_base
    print 'p0_std', p0_std
    p0 = np.random.normal(p0_base, p0_std, size=[ntemps, nwalkers, ndim])
    print 'p0', p0
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
    if NUFIT:
        outfile += '_nufit'
    if FLAT:
	outfile += '_flat'

    np.save(outfile+'.npy', samples)

    print "Making triangle plots"
    chainer_plot.plot(
        infile=outfile+'.npy',
        angles=True,
        nufit=NUFIT,
        outfile=outfile+'_angles.pdf',
        bestfit_ratio=BESTFIT,
        sigma_ratio=SIGMA
    )
    chainer_plot.plot(
        infile=outfile+'.npy',
        angles=False,
        nufit=NUFIT,
        outfile=outfile+'.pdf',
        bestfit_ratio=BESTFIT,
        sigma_ratio=SIGMA
    )
    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
