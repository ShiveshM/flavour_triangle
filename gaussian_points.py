#! /usr/bin/env python
"""
Generate a sample of flavour ratio points with gaussian fluctuations
"""

from __future__ import absolute_import, division

import argparse

import numpy as np

from pisa.utils.log import logging, set_verbosity


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--ratio', type=int, nargs=3, default=[1, 1, 1],
        help='Set the flavour ratio'
    )
    parser.add_argument(
        '--sigma', type=float, default=0.2,
        help='1 sigma deviation'
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

    np.random.seed(args.seed)

    with open(args.outfile, 'w') as f:
        logging.debug('Writing to file {0}'.format(args.outfile))
        for i in xrange(args.number):
            logging.trace('Point #{0}'.format(i))
            mean = args.ratio[:-1]
            sigma = np.identity(len(mean)) * args.sigma
            logging.trace('sigma = \n{0}'.format(sigma))

            sucess = False
            while not sucess:
                points = np.random.multivariate_normal(
                    mean=mean, cov=sigma
                )
                if np.any(points < 0) or np.sum(points) > 1:
                    logging.trace('Variation {0} out of range, '
                                  'retrying'.format(points))
                    continue
                sucess = True

            flav_comp = points.tolist() + [1 - np.sum(points)]
            for x in flav_comp:
                f.write('{0:.6f} '.format(x))
            f.write('\n')


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
