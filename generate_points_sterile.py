#! /usr/bin/env python
"""
Generate a sample of flavour ratio points with fluctuations from a
unitary 4x4 mixing matrix.
"""

from __future__ import absolute_import, division

import argparse

import numpy as np

from pisa.core.prior import Prior
from pisa.core.param import Param, ParamSet
from pisa.utils.log import logging, set_verbosity
from pisa.utils.profiler import profile


def compute_flavour_ratio(initial_ratio, mixing_matrix):
    """Compute the observed flavour ratio assuming decoherence."""
    logging.debug('Entering compute_flavour_ratio')
    composition = np.einsum(
        'ai, bi, a -> b', mixing_matrix**2, mixing_matrix**2, initial_ratio
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


def randomise_params(mm_params, non_unitarity, maxtrials):
    """Randomise a set of mixing matrix values according to a gaussian
    prior.
    """
    logging.debug('Entering randomise_params')
    success = False
    i = 0
    while not success:
        if i == maxtrials:
            raise ValueError('Maxtrials has been reached')
        randomise_param(mm_params[0])
        randomise_param(mm_params[1])
        i +=1
        mm2_value_2 = non_unitarity - (mm_params[0].value.m**2 +
                                       mm_params[1].value.m**2)
        if mm2_value_2 < 0.:
            logging.trace('Variation {0} out of range for param {1}, '
                          'retrying'.format(mm2_value_2, mm_params[2].name))
            continue
        mm2_value = np.sqrt(mm2_value_2)
        try:
            mm_params[2].value = mm2_value
        except ValueError:
            # logging.debug('0 = {0} , 1 = {1} , nu = {2} , mm2_value = {3}'.format(
            #     mm_params[0].value.m, mm_params[1].value.m, non_unitarity, mm2_value
            # ))
            logging.trace('Variation {0} out of range for param {1}, '
                          'retrying'.format(mm2_value, mm_params[2].name))
            continue
        success = True


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
            continue
        try:
            for row in mm_params:
                randomise_params(row, non_unitarity, maxtrials)
        except ValueError:
            logging.debug('Maxtrials {0} reached for non-unitarity = {1}, '
                          'trying another value of non-unitarity'.format(
                              maxtrials, non_unitarity
                          ))
            continue
        success = True
    logging.debug('randomised mixing_matrix_paramset '
                  '=\n{0}'.format(mixing_matrix_paramset))


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
            randomise_paramset(
                mixing_matrix_paramset=mm_paramset,
                central_nonunitarity=args.central_nonunitarity,
                sigma_nonunitarity=args.sigma_nonunitarity,
                maxtrials=args.maxtrials
            )

            mm_matrix = create_matrix(
                mixing_matrix_paramset=mm_paramset
            )

            flav_comp = compute_flavour_ratio(
                args.initial_ratio, mm_matrix
            )

            for x in flav_comp:
                f.write('{0:.6f} '.format(x))
            f.write('{0:.6f}'.format(np.sum(flav_comp)))
            f.write('\n')


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
