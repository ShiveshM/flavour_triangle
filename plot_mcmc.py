#! /usr/bin/env python
"""
Plot chains from MCMC
"""

from __future__ import absolute_import, division

import argparse

import numpy as np
import matplotlib
import corner


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--chains', type=str, required=True,
        help='Path to chains'
    )
    parser.add_argument(
        '--outfile', type=str, default='./untitled.npy',
        help='Path to output chains'
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    print 'Making plot'
    samples = np.load(args.chains)
    fig = corner.corner(samples)
    fig.savefig(args.outfile)
    print 'DONE!'


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
