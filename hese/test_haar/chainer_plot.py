#! /usr/bin/env python
"""
From an MCMC chains file, make a triangle plot.
"""

from __future__ import absolute_import, division

import sys
sys.path.append('/users/mandalia/Documents/flavour_triangle/hese/test_haar')

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import rc, rcParams

import getdist
from getdist import plots
from getdist import mcsamples

import mcmc_scan


rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})
cols = ['#29A2C6','#FF6D31','#FFCB18','#73B66B','#EF597B', '#333333']

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 18}


def plot(infile, angles, outfile, bestfit_ratio=None, sigma_ratio=None):
    """Make the triangle plot"""
    if not angles:
        labels=[r'\phi_e', r'\phi_\mu', r'\phi_\tau']
    else:
        labels=[r'sin^4(\theta)', r'cos(2\phi)']
    print labels

    if not angles:
        ranges = []
        for x in xrange(len(labels)):
            ranges.append([0, 1])
    else:
        ranges = [(0., 1.), (-1., 1.)]

    raw = np.load(infile)
    if not angles:
        Tchain = np.array(map(mcmc_scan.coord_transform, raw))
    else:
        Tchain = raw

    if bestfit_ratio is not None and sigma_ratio is not None:
        label = 'Bestfit ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nSigma = {3:.3f}'.format(
            bestfit_ratio[0], bestfit_ratio[1], bestfit_ratio[2], sigma_ratio
        )
    else:
        label = None

    Tsample = mcsamples.MCSamples(
        samples=Tchain, labels=labels, ranges=ranges
    )

    Tsample.updateSettings({'contours': [0.90, 0.99]})
    Tsample.num_bins_2D=500
    Tsample.fine_bins_2D=500
    Tsample.smooth_scale_2D=0.03

    g = plots.getSubplotPlotter()
    g.settings.num_plot_contours = 2
    g.settings.axes_fontsize = 10
    g.settings.figure_legend_frame = False
    g.triangle_plot(
        [Tsample], filled=True,
    )
    mpl.pyplot.figtext(0.6, 0.7, label, fontsize=15)
    print 'outfile = {0}'.format(outfile)
    g.export(outfile)


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--infile', type=str, required=True,
        help='Path to MCMC chains'
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
        '--bestfit-ratio', type=int, nargs=3, required=False,
        help='Set the bestfit flavour ratio'
    )
    parser.add_argument(
        '--sigma-ratio', type=float, required=False,
        help='Set the 1 sigma for the flavour ratio'
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
