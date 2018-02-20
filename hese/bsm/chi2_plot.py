#! /usr/bin/env python
"""
From a sampled space of likelihood values, make a chi2 plot.
"""

from __future__ import absolute_import, division

import os, sys
import errno
import argparse

import numpy as np
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 1000
mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import rc, rcParams

import h5py

rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})
cols = ['#29A2C6','#FF6D31','#FFCB18','#73B66B','#EF597B', '#333333']

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 18}


def plot(infile, angles, outfile, measured_ratio, sigma_ratio, fix_sfr,
         fix_mixing, fix_scale, source_ratio, scale, dimension, energy,
         scale_bounds):
    """Make the chi2 plot."""
    raw = np.array(h5py.File(infile)['df']).T
    print 'raw.shape', raw.shape

    scale = raw[-2]
    argmin_scale = np.argmin(scale)

    llh = raw[-1]
    null_llh = llh[argmin_scale]
    print 'null LLH', raw.T[argmin_scale]

    chi2 = -2 * (llh - null_llh)
    print 'chi2', chi2

    x_label = r'${\rm log}_{10}\Lambda$'
    y_label = r'$-2\Delta LLH$'
    if fix_sfr:
        if fix_scale:
            plot_label = 'Source flavour ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nIC observed flavour ratio = [{3:.2f}, {4:.2f}, {5:.2f}]\nSigma = {6:.3f}\nDimension = {7}\nEnergy = {8} GeV\nScale = {9}'.format(
                source_ratio[0], source_ratio[1], source_ratio[2],
                measured_ratio[0], measured_ratio[1], measured_ratio[2], sigma_ratio,
                dimension, int(energy), scale
            )
        else:
            plot_label = 'Source flavour ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nIC observed flavour ratio = [{3:.2f}, {4:.2f}, {5:.2f}]\nSigma = {6:.3f}\nDimension = {7}\nEnergy = {8} GeV'.format(
                source_ratio[0], source_ratio[1], source_ratio[2],
                measured_ratio[0], measured_ratio[1], measured_ratio[2], sigma_ratio,
                dimension, int(energy)
            )
    else:
        if fix_scale:
	    plot_label = 'IC observed flavour ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nSigma = {3:.3f}\nDimension = {4}\nEnergy = {5} GeV\nScale = {6}'.format(
		measured_ratio[0], measured_ratio[1], measured_ratio[2], sigma_ratio,
		dimension, int(energy), scale
	    )
	else:
	    plot_label = 'IC observed flavour ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nSigma = {3:.3f}\nDimension = {4}\nEnergy = {5} GeV'.format(
		measured_ratio[0], measured_ratio[1], measured_ratio[2], sigma_ratio,
		dimension, int(energy)
	    )

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)

    ax.tick_params(axis='x', labelsize=11)
    ax.tick_params(axis='y', labelsize=11)
    ax.set_xlim(np.log10(scale_bounds))

    fig.text(0.5, -0.01, x_label, ha='center', fontsize=18)
    fig.text(0.01, 0.5, y_label, va='center', rotation='vertical', fontsize=18)
    fig.text(0.2, 0.7, plot_label, fontsize=8)

    ax.scatter(scale, chi2, s=1, c='b', marker='.')

    print 'outfile = {0}'.format(outfile)
    try:
        os.makedirs(outfile[:-len(os.path.basename(outfile))])
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(outfile[:-len(os.path.basename(outfile))]):
            pass
        else:
            raise
    fig.savefig(outfile, bbox_inches='tight', dpi=150)


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--infile', type=str, required=True,
        help='Path to liklihood'
    )
    parser.add_argument(
        '--angles', default=False, action='store_true',
        help='Plot in terms of mixing angles'
    )
    parser.add_argument(
        '--outfile', type=str, default='./untitled.pdf',
        help='Path to output plot'
    )
    parser.add_argument(
        '--measured-ratio', type=int, nargs=3, required=False,
        help='Set the measured flavour ratio'
    )
    parser.add_argument(
        '--sigma-ratio', type=float, required=False,
        help='Set the 1 sigma for the flavour ratio'
    )
    parser.add_argument(
        '--fix-sfr', action='store_true',
        help='Fix the source flavour ratio'
    )
    parser.add_argument(
        '--fix-mixing', action='store_true',
        help='Fix the new physics mixing values to a single term, s_12^2'
    )
    parser.add_argument(
        '--fix-scale', action='store_true',
        help='Fix the new physics scale'
    )
    parser.add_argument(
        '--source-ratio', type=int, nargs=3, default=[2, 1, 0],
        help='Set the source flavour ratio for the case when you want to fix it'
    )
    parser.add_argument(
        '--scale', type=float, required=False,
        help='Fix the scale to this value'
    )
    parser.add_argument(
        '--dimension', type=int, default=3, help='Dimension'
    )
    parser.add_argument(
        '--energy', type=float, default=1000, help='Energy'
    )
    parser.add_argument(
        '--scale-bounds', type=float, nargs=2,
        help='Upper and lower limits to plot the new physics scale'
    )
    args = parser.parse_args()
    return args


def main():
    args = vars(parse_args())
    plot(**args)

    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
