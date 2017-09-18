#! /usr/bin/env python
"""
Generate a sample of flavour ratio points with fluctuations from a non-unitary
mixing matrix with the contraint that the matrix is a subset of a unitary 4x4
matrix.
"""

from __future__ import absolute_import, division

import argparse

import numpy as np
import scipy.linalg

from pisa.core.prior import Prior
from pisa.core.param import Param, ParamSet
from pisa.utils.log import logging, set_verbosity
from pisa.utils.profiler import profile


def get_std_components(mm_paramset):
    """Return the standard mixing matrix parameters as a 2D numpy array."""
    logging.trace('Entering get_std_components')
    Ue1 = mm_paramset.Ue1
    Ue2 = mm_paramset.Ue2
    Ue3 = mm_paramset.Ue3
    Um1 = mm_paramset.Um1
    Um2 = mm_paramset.Um2
    Um3 = mm_paramset.Um3
    Ut1 = mm_paramset.Ut1
    Ut2 = mm_paramset.Ut2
    Ut3 = mm_paramset.Ut3
    return np.array([(Ue1, Ue2, Ue3), (Um1, Um2, Um3), (Ut1, Ut2, Ut3)])


def get_all_components(mm_paramset):
    """Return the mixing matrix parameters as a 2D numpy array."""
    logging.trace('Entering get_all_components')
    Ue1 = mm_paramset.Ue1
    Ue2 = mm_paramset.Ue2
    Ue3 = mm_paramset.Ue3
    Ue4 = mm_paramset.Ue4
    Um1 = mm_paramset.Um1
    Um2 = mm_paramset.Um2
    Um3 = mm_paramset.Um3
    Um4 = mm_paramset.Um4
    Ut1 = mm_paramset.Ut1
    Ut2 = mm_paramset.Ut2
    Ut3 = mm_paramset.Ut3
    Ut4 = mm_paramset.Ut4
    Us1 = mm_paramset.Us1
    Us2 = mm_paramset.Us2
    Us3 = mm_paramset.Us3
    Us4 = mm_paramset.Us4
    return np.array([(Ue1, Ue2, Ue3, Ue4), (Um1, Um2, Um3, Um4),
                     (Ut1, Ut2, Ut3, Ut4), (Us1, Us2, Us3, Us4)])


def compute_flavour_ratio(initial_ratio, mixing_matrix):
    """Compute the observed flavour ratio assuming decoherence."""
    logging.debug('Entering compute_flavour_ratio')
    composition = np.einsum(
        'ai, bi, a -> b', abs(mixing_matrix)**2, abs(mixing_matrix)**2, initial_ratio
    )
    ratio = composition / np.sum(initial_ratio)
    logging.trace('Computed flavour ratio = {0}'.format(ratio))
    return ratio


def test_nuni(x):
    """Test the non-unitarity of a matrix."""
    logging.trace('Entering test_uni')
    logging.debug('Unitarity test:\n{0}'.format(abs(np.dot(x, x.conj().T))))


def create_matrix(mm_paramset):
    """Create a matrix from a ParamSet of mixing matrices."""
    logging.debug('Entering create_matrix')
    def get_m(a):
        return a.value
    v_get_m = np.vectorize(get_m)
    matrix = v_get_m(get_all_components(mm_paramset)).astype(np.complex128)
    logging.trace('matrix = \n{0}'.format(matrix))
    return matrix


def make_unitary(x):
    """Create a unitary matrix from a given matrix."""
    q, r = scipy.linalg.qr(x)
    d = r.diagonal()
    q *= d/abs(d)
    return q


def set_nonstd_params(mm_matrix, nuni):
    """Set the nonstd params given a non-unitarity value"""
    logging.debug('Entering set_nonstd_params')
    def maths(z):
        z[-1] = np.sqrt(1. + 0j - nuni - np.sum(
            map(lambda x: abs(x)**2, z[:-1])
        ))

    map(maths, mm_matrix[:-1])
    map(maths, mm_matrix.T[:-1])
    # TODO(shivesh): check if requiring the last element to add to nuni both
    # row and columns wise matters
    maths(mm_matrix[-1])
    logging.debug('Set nonstd_params, mm =\n{0}'.format(mm_matrix))
    logging.debug('abs =\n{0}'.format(abs(mm_matrix)))


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


def initialise_paramset():
    """Initialise mixing matrix with priors from arxiv 1508.05095."""
    logging.debug('Entering initialise_paramset')
    # TODO(shivesh): these are obtained by eye!

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

    mm_pm = []
    mm_pm.append(gauss(name=r'Ue1', mean=0.82, stddev=0.02))
    mm_pm.append(gauss(name=r'Ue2', mean=0.55, stddev=0.02))
    mm_pm.append(gauss(name=r'Ue3', mean=0.15, stddev=0.01))
    mm_pm.append(gauss(name=r'Ue4', mean=0.00, stddev=0.50))
    mm_pm.append(gauss(name=r'Um1', mean=0.38, stddev=0.13))
    mm_pm.append(gauss(name=r'Um2', mean=0.54, stddev=0.10))
    mm_pm.append(gauss(name=r'Um3', mean=0.75, stddev=0.10))
    mm_pm.append(gauss(name=r'Um4', mean=0.00, stddev=0.50))
    mm_pm.append(gauss(name=r'Ut1', mean=0.41, stddev=0.14))
    mm_pm.append(gauss(name=r'Ut2', mean=0.63, stddev=0.13))
    mm_pm.append(gauss(name=r'Ut3', mean=0.64, stddev=0.11))
    mm_pm.append(gauss(name=r'Ut4', mean=0.00, stddev=0.50))
    mm_pm.append(gauss(name=r'Us1', mean=0.00, stddev=0.50))
    mm_pm.append(gauss(name=r'Us2', mean=0.00, stddev=0.50))
    mm_pm.append(gauss(name=r'Us3', mean=0.00, stddev=0.50))
    mm_pm.append(gauss(name=r'Us4', mean=0.00, stddev=0.50))

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
        '--nuni', type=float, default=0.1,
        help='Amount of non-unitarity to impose'
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

    bufsize = 0
    with open(args.outfile, 'w', bufsize) as f:
        logging.debug('Writing to file {0}'.format(args.outfile))
        for i in xrange(args.number):
            logging.trace('Point #{0}'.format(i))

            std_mm_params = get_std_components(mm_paramset)
            success = False
            while not success:
                for y in std_mm_params:
                    for x in y:
                        randomise_param(x)

                mm_matrix = create_matrix(
                    mm_paramset=mm_paramset
                )

                set_nonstd_params(
                    mm_matrix=mm_matrix,
                    nuni=args.nuni
                )

                u_mm_matrix = make_unitary(mm_matrix)
                if args.v > 1:
                    test_nuni(u_mm_matrix)

                fail = False
                for i in np.ndindex(std_mm_params.shape):
                    try:
                        std_mm_params[i].value = abs(u_mm_matrix[i])
                    except ValueError:
                        logging.trace(
                            'Unitary matrix element {0} has value {1:.2f} '
                            'which is out of bounds, retrying'.format(
                                std_mm_params[i].name, abs(u_mm_matrix[i])
                            )
                        )
                        fail = True
                        break
                if fail: continue
                success = True

            flav_comp = compute_flavour_ratio(
                args.initial_ratio + [0], u_mm_matrix
            )

            three_flav_comp = flav_comp[:-1]

            for x in three_flav_comp:
                f.write('{0:.6f} '.format(x))
            f.write('\n')

main.__doc__ = __doc__


if __name__ == '__main__':
    main()
