#! /usr/bin/env python
"""
Generate a sample of flavour ratio points with fluctuations from a non-unitary
mixing matrix.
"""

from __future__ import absolute_import, division

import argparse

import numpy as np
from pisa import ureg, Q_
from pisa.core.prior import Prior
from pisa.core.param import Param, ParamSet
from pisa.utils.log import logging, set_verbosity


def compute_flavour_ratio(initial_ratio, mixing_matrix):
    """Compute the observed flavour ratio assuming decoherence."""
    logging.debug('Entering compute_flavour_ratio')
    composition = np.einsum(
        'ai, bi, a -> b', mixing_matrix**2, mixing_matrix**2, initial_ratio
    )
    ratio = composition / np.sum(initial_ratio)
    logging.trace('Computed flavour ratio = {0}'.format(ratio))
    return ratio


def randomise_paramset(mixing_matrix_paramset):
    """Randomise the mixing matrix values according to a gaussian prior."""
    logging.debug('Entering randomise_paramset')
    # TODO(shivesh): gaussian only
    logging.trace('initial mixing_matrix_paramset '
                  '=\n{0}'.format(mixing_matrix_paramset))
    for param in mixing_matrix_paramset:
        sucess = False
        while not sucess:
            rndm = np.random.normal(param.nominal_value, param.prior.stddev)
            try:
                param.validate_value(rndm)
            except ValueError:
                continue
            sucess = True
        param.value = rndm
    logging.trace('randomised mixing_matrix_paramset '
                  '=\n{0}'.format(mixing_matrix_paramset))


def create_matrix(mixing_matrix_paramset):
    """Create a matrix from a ParamSet of mixing matrices."""
    logging.debug('Entering create_matrix')
    Ue1 = mixing_matrix_paramset.Ue1.value.m
    Ue2 = mixing_matrix_paramset.Ue2.value.m
    Ue3 = mixing_matrix_paramset.Ue3.value.m
    Um1 = mixing_matrix_paramset.Um1.value.m
    Um2 = mixing_matrix_paramset.Um2.value.m
    Um3 = mixing_matrix_paramset.Um3.value.m
    Ut1 = mixing_matrix_paramset.Ut1.value.m
    Ut2 = mixing_matrix_paramset.Ut2.value.m
    Ut3 = mixing_matrix_paramset.Ut3.value.m
    matrix = np.array([[Ue1, Ue2, Ue3], [Um1, Um2, Um3], [Ut1, Ut2, Ut3]])
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
        '--number', type=int, default=100,
        help='Number of points to sample'
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
                mixing_matrix_paramset=mm_paramset
            )

            mm_matrix = create_matrix(
                mixing_matrix_paramset=mm_paramset
            )

            flav_comp = compute_flavour_ratio(
                args.initial_ratio, mm_matrix
            )

            for x in flav_comp:
                f.write('{0:.6f} '.format(x))
            f.write('\n')

main.__doc__ = __doc__


if __name__ == '__main__':
    main()
