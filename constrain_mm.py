#! /usr/bin/env python
"""
From an dataset of flavour traingle points, find the exclusions for the
mixing matrix elements.
"""

from __future__ import absolute_import, division

import argparse

import numpy as np

from pisa.core.prior import Prior
from pisa.core.param import Param, ParamSet
from pisa.utils.log import logging, set_verbosity
from pisa.utils.profiler import profile


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--infile', type=str, required=True,
        help='Input file with flavour ratios'
    )
    parser.add_argument(
        '--bins', type=int, default=5,
        help='Number of bins'
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

    logging.info('Loading data from {0}'.format(args.infile))
    flavour_list = np.genfromtxt(args.infile)
    logging.debug('Finished loading data from {0}'.format(args.infile))
    logging.debug('Data shape = {0}'.format(flavour_list.shape))

    binning = np.linspace(0, 1, args.bins + 1)
    hist, _ = np.histogramdd(
        sample=flavour_list,
        bins=[binning] * flavour_list.shape[1],
    )
    logging.trace('Hist shape = {0}'.format(hist.shape))


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
