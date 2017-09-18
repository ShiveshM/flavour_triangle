#! /usr/bin/env python
"""
Generate a sample of flavour ratio points with fluctuations from a non-unitary
mixing matrix.
"""

from __future__ import absolute_import, division

import argparse

import numpy as np
import scipy.linalg

from pisa.core.prior import Prior
from pisa.core.param import Param, ParamSet
from pisa.utils.log import logging, set_verbosity
from pisa.utils.profiler import profile


def compute_flavour_ratio(initial_ratio, mixing_matrix):
    """Compute the observed flavour ratio assuming decoherence."""
    logging.debug('Entering compute_flavour_ratio')
    composition = np.einsum(
        'ai, bi, a -> b', abs(mixing_matrix)**2, abs(mixing_matrix)**2, initial_ratio
    )
    ratio = composition / np.sum(initial_ratio)
    logging.trace('Computed flavour ratio = {0}'.format(ratio))
    return ratio


def get_mm_params(mixing_matrix_paramset):
    """Return the mixing matrix parameters as a 2D tuple."""
    logging.trace('Entering get_mm_params')
    Ue1 = mixing_matrix_paramset.Ue1
    Ue2 = mixing_matrix_paramset.Ue2
    Ue3 = mixing_matrix_paramset.Ue3
    Um1 = mixing_matrix_paramset.Um1
    Um2 = mixing_matrix_paramset.Um2
    Um3 = mixing_matrix_paramset.Um3
    Ut1 = mixing_matrix_paramset.Ut1
    Ut2 = mixing_matrix_paramset.Ut2
    Ut3 = mixing_matrix_paramset.Ut3
    return ((Ue1, Ue2, Ue3), (Um1, Um2, Um3), (Ut1, Ut2, Ut3))


def randomise_param(mm_param):
    """Randomise a single mixing matrix value according to its bounds."""
    logging.trace('Entering randomise_param')
    success = False
    while not success:
        rndm = np.random.normal(mm_param.nominal_value, mm_param.prior.stddev)
        try:
            mm_param.value = rndm
        except ValueError:
            logging.trace('Variation {0} out of range for param {1}, '
                          'retrying'.format(rndm, mm_param.name))
            continue
        success = True


def make_unitary(x):
    """Create a unitary matrix from a given matrix."""
    q, r = scipy.linalg.qr(x)
    d = r.diagonal()
    q *= d/abs(d)
    return q


def gauss(name, mean, stddev, **kwargs):
    logging.trace('Entering gauss: name = {0}, mean = {1}, stddev = '
                  '{1}'.format(name, mean, stddev))
    prior= Prior(kind='gaussian', mean=mean, stddev=stddev)
    range = 3 * np.array([-stddev, stddev]) + mean
    if range[0] < 0: range[0] = 0.01
    if range[1] > 1: range[1] = 0.99
    return Param(
        name=name, value=mean, prior=prior, range=range, is_fixed=False,
        is_discrete=False, **kwargs
    )


def randomise_params(mm_params, non_unitarity, maxtrials):
    """Randomise a set of mixing matrix values according to a gaussian
    prior. Constraints are given by having both the rows and the columns being
    summed to a specific value of non-unitarity and also by the 3 sigma C.L.
    given by Parke's paper.
    """
    logging.debug('Entering randomise_params')
    success = False
    i = 0
    while not success:
        if i == maxtrials:
            raise ValueError('Maxtrials has been reached')
        for x in mm_params:
            for y in x:
                randomise_param(y)
        i +=1
        z = []
        for y in mm_params:
            z += [x.value.m for x in y]
            z.append(np.sqrt(non_unitarity + 0*1j - np.sum(map(lambda x: x.value.m**2, y))))
        for y in np.array(mm_params).T:
            z.append(np.sqrt(non_unitarity + 0*1j - np.sum(map(lambda x: x.value.m**2, y))))
        z.append(np.sqrt(non_unitarity - (z[12]**2 + z[13]**2 + z[14]**2)))
        z = np.array(z).reshape(4, 4)

        z_u = make_unitary(z)
        mm_params_4d = np.zeros((4, 4))
        try:
            for x in xrange(4):
                for y in xrange(4):
                    if x < 3 and y < 3:
                        mm_params[x][y].value = abs(z_u[x][y])
                    mm_params_4d[x][y] = z_u[x][y]
        except ValueError:
            # logging.debug('0 = {0} , 1 = {1} , nu = {2} , mm2_value = {3}'.format(
            #     mm_params[0].value.m, mm_params[1].value.m, non_unitarity, mm2_value
            # ))
            logging.trace('Variation {0} out of range'.format(mm_params_4d))
            continue
        success = True

    return mm_params_4d


@profile
def randomise_paramset(mixing_matrix_paramset, central_nonunitarity,
                       sigma_nonunitarity, maxtrials):
    """Randomise all mixing matrix values according to a gaussian prior."""
    logging.debug('Entering randomise_paramset')
    # TODO(shivesh): gaussian only
    logging.trace('initial mixing_matrix_paramset '
                  '=\n{0}'.format(mixing_matrix_paramset))
    mm_params = get_mm_params(mixing_matrix_paramset)
    success = False
    while not success:
        non_unitarity = np.random.normal(
            central_nonunitarity, sigma_nonunitarity
        )
        logging.trace('Setting non-unitarity to {0}'.format(non_unitarity))
        if non_unitarity < 0.:
            logging.trace('Non unitarity set to negative value, trying '
                          'another value of non-unitarity.')
            continue
        try:
            mm_params_4d = randomise_params(mm_params, non_unitarity,
                                            maxtrials)
        except ValueError:
            logging.debug('Maxtrials {0} reached for non-unitarity = {1}, '
                          'trying another value of non-unitarity'.format(
                              maxtrials, non_unitarity
                          ))
            continue
        success = True

    s_mixing_matrix_paramset = []
    for x in mixing_matrix_paramset:
        s_mixing_matrix_paramset.append(x)
    s_mixing_matrix_paramset.append(gauss(name=r'Ue4', mean=mm_params_4d[0][3],
                                          stddev=0.01))
    s_mixing_matrix_paramset.append(gauss(name=r'Um4', mean=mm_params_4d[1][3],
                                          stddev=0.01))
    s_mixing_matrix_paramset.append(gauss(name=r'Ut4', mean=mm_params_4d[2][3],
                                          stddev=0.01))
    s_mixing_matrix_paramset.append(gauss(name=r'Us1', mean=mm_params_4d[3][0],
                                          stddev=0.01))
    s_mixing_matrix_paramset.append(gauss(name=r'Us2', mean=mm_params_4d[3][1],
                                          stddev=0.01))
    s_mixing_matrix_paramset.append(gauss(name=r'Us3', mean=mm_params_4d[3][2],
                                          stddev=0.01))
    s_mixing_matrix_paramset.append(gauss(name=r'Us4', mean=mm_params_4d[3][3],
                                          stddev=0.01))
    s_mixing_matrix_paramset = ParamSet(s_mixing_matrix_paramset)

    logging.debug('randomised mixing_matrix_paramset '
                  '=\n{0}'.format(s_mixing_matrix_paramset))
    return s_mixing_matrix_paramset


def create_matrix(mixing_matrix_paramset):
    """Create a matrix from a ParamSet of mixing matrices."""
    logging.debug('Entering create_matrix')
    def get_m(a):
        return a.value.m
    v_get_m = np.vectorize(get_m)
    mm_params = v_get_m(get_mm_params(mixing_matrix_paramset))
    matrix = np.array(mm_params)
    logging.trace('matrix = \n{0}'.format(matrix))
    return matrix


def initialise_paramset():
    """Initialise mixing matrix with priors from arxiv 1508.05095."""
    logging.debug('Entering initialise_paramset')
    # TODO(shivesh): use splines instead of gaussian

    # TODO(shivesh): these are obtained by eye!
    mm_pm = []
    mm_pm.append(gauss(name=r'Ue1', mean=0.82, stddev=0.02))
    mm_pm.append(gauss(name=r'Ue2', mean=0.55, stddev=0.02))
    mm_pm.append(gauss(name=r'Ue3', mean=0.15, stddev=0.01))
    mm_pm.append(gauss(name=r'Um1', mean=0.38, stddev=0.13))
    mm_pm.append(gauss(name=r'Um2', mean=0.54, stddev=0.10))
    mm_pm.append(gauss(name=r'Um3', mean=0.75, stddev=0.10))
    mm_pm.append(gauss(name=r'Ut1', mean=0.41, stddev=0.14))
    mm_pm.append(gauss(name=r'Ut2', mean=0.63, stddev=0.13))
    mm_pm.append(gauss(name=r'Ut3', mean=0.64, stddev=0.11))

    for pm in mm_pm:
        logging.trace('pm = {0}'.format(pm))
    return ParamSet(mm_pm)


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--initial-ratio', type=int, nargs=3, default=[1, 2, 0],
        help='Set the initial flavour ratio at production'
    )
    parser.add_argument(
        '--outfile', type=str, default='./untitled.txt',
        help='Path to output file'
    )
    parser.add_argument(
        '--central-nonunitarity', type=float, default=1.0,
        help='central value for non-unitarity'
    )
    parser.add_argument(
        '--sigma-nonunitarity', type=float, default=0.2,
        help='1 sigma limit for non-unitarity'
    )
    parser.add_argument(
        '--number', type=int, default=100,
        help='Number of points to sample'
    )
    parser.add_argument(
        '--maxtrials', type=int, default=100,
        help='Maximum number of trials to try per non-unitarity value'
    )
    parser.add_argument(
        '--seed', type=int, default=99,
        help='Set the random seed value'
    )
    parser.add_argument(
        '-v', action='count',
        help='Set verbosity level; repeat -v for higher level.'
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    set_verbosity(args.v)

    mm_paramset = initialise_paramset()

    np.random.seed(args.seed)

    with open(args.outfile, 'w') as f:
        logging.debug('Writing to file {0}'.format(args.outfile))
        for i in xrange(args.number):
            logging.trace('Point #{0}'.format(i))
            s_mm_paramset = randomise_paramset(
                mixing_matrix_paramset=mm_paramset,
                central_nonunitarity=args.central_nonunitarity,
                sigma_nonunitarity=args.sigma_nonunitarity,
                maxtrials=args.maxtrials
            )

            print 'HERE', s_mm_paramset
            mm_matrix = create_matrix(
                mixing_matrix_paramset=s_mm_paramset
            )

            flav_comp = compute_flavour_ratio(
                args.initial_ratio + [0], mm_matrix
            )

            for x in flav_comp:
                f.write('{0:.6f} '.format(x))
            f.write('{0:.6f}'.format(np.sum(flav_comp)))
            f.write('\n')


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
