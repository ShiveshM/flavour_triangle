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
from matplotlib.ticker import FormatStrFormatter

import h5py
import pandas as pd
from pandas.tools.plotting import scatter_matrix

rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})
cols = ['#29A2C6','#FF6D31','#FFCB18','#73B66B','#EF597B', '#333333']

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 18}


CRIT_90_CHI2 = {
    2 : 4.605,
    3 : 6.251,
    4 : 7.779,
    5 : 9.236,
    6 : 10.645,
    7 : 12.017,
    8 : 13.362
}


def plot(infile, angles, outfile, measured_ratio, sigma_ratio, fix_sfr,
         fix_mixing, fix_scale, source_ratio, scale, dimension, energy,
         scale_bounds):
    """Make the limit plot."""
    raw = np.array(h5py.File(infile)['df']).T
    print 'raw.shape', raw.shape

    if not angles:
        if fix_mixing:
            labels = [r'{\rm log}_{10}\Lambda', r'\phi_e', r'\phi_\mu', r'\phi_\tau']
        elif fix_sfr:
            if fix_scale:
                labels = [r'\mid \tilde{U}_{e1} \mid', r'\mid \tilde{U}_{e2} \mid', r'\mid \tilde{U}_{e3} \mid', \
                          r'\mid \tilde{U}_{\mu1} \mid', r'\mid \tilde{U}_{\mu2} \mid', r'\mid \tilde{U}_{\mu3} \mid', \
                          r'\mid \tilde{U}_{\tau1} \mid', r'\mid \tilde{U}_{\tau2} \mid', r'\mid \tilde{U}_{\tau3} \mid']
            else:
                labels = [r'\mid \tilde{U}_{e1} \mid', r'\mid \tilde{U}_{e2} \mid', r'\mid \tilde{U}_{e3} \mid', \
                          r'\mid \tilde{U}_{\mu1} \mid', r'\mid \tilde{U}_{\mu2} \mid', r'\mid \tilde{U}_{\mu3} \mid', \
                          r'\mid \tilde{U}_{\tau1} \mid', r'\mid \tilde{U}_{\tau2} \mid', r'\mid \tilde{U}_{\tau3} \mid', \
                          r'{\rm log}_{10}(\Lambda)']
        else:
            if fix_scale:
                labels = [r'\mid \tilde{U}_{e1} \mid', r'\mid \tilde{U}_{e2} \mid', r'\mid \tilde{U}_{e3} \mid', \
                          r'\mid \tilde{U}_{\mu1} \mid', r'\mid \tilde{U}_{\mu2} \mid', r'\mid \tilde{U}_{\mu3} \mid', \
                          r'\mid \tilde{U}_{\tau1} \mid', r'\mid \tilde{U}_{\tau2} \mid', r'\mid \tilde{U}_{\tau3} \mid', \
                          r'{\rm log}_{10}\Lambda', r'\phi_e', r'\phi_\mu', r'\phi_\tau']
            else:
                labels = [r'\mid \tilde{U}_{e1} \mid', r'\mid \tilde{U}_{e2} \mid', r'\mid \tilde{U}_{e3} \mid', \
                          r'\mid \tilde{U}_{\mu1} \mid', r'\mid \tilde{U}_{\mu2} \mid', r'\mid \tilde{U}_{\mu3} \mid', \
                          r'\mid \tilde{U}_{\tau1} \mid', r'\mid \tilde{U}_{\tau2} \mid', r'\mid \tilde{U}_{\tau3} \mid', \
                          r'\phi_e', r'\phi_\mu', r'\phi_\tau']
    else:
        if fix_sfr:
            if fix_mixing:
                assert 0
                labels=[r'\tilde{s}_{12}^2', r'{\rm log}_{10}\Lambda']
            elif fix_scale:
                labels=[r'\tilde{s}_{12}^2', r'\tilde{c}_{13}^4',
                        r'\tilde{s}_{23}^2', r'\tilde{\delta_{CP}}']
            else:
                labels=[r'\tilde{s}_{12}^2', r'\tilde{c}_{13}^4',
                        r'\tilde{s}_{23}^2', r'\tilde{\delta_{CP}}',
                        r'{\rm log}_{10}\Lambda']
        else:
            if fix_mixing:
                labels=[r'{\rm log}_{10}\Lambda', r'sin^4(\phi)', r'cos(2\psi)']
            elif fix_scale:
                labels=[r'\tilde{s}_{12}^2', r'\tilde{c}_{13}^4',
                        r'\tilde{s}_{23}^2', r'\tilde{\delta_{CP}}',
                        r'sin^4(\phi)', r'cos(2\psi)']
            else:
                labels=[r'\tilde{s}_{12}^2', r'\tilde{c}_{13}^4',
                        r'\tilde{s}_{23}^2', r'\tilde{\delta_{CP}}',
                        r'{\rm log}_{10}\Lambda', r'sin^4(\phi)', r'cos(2\psi)']
    labels = ['$'+x+'$' for x in labels]
    print 'labels', labels

    scale = raw[-2]
    argmin_scale = np.argmin(scale)

    llh = raw[-1]
    null_llh = llh[argmin_scale]
    print 'null LLH', raw.T[argmin_scale]

    chi2 = -2 * (llh - null_llh)
    print 'chi2', chi2
    print 'critical 90% chi2 value', CRIT_90_CHI2[len(labels)]

    chi2_mask = (chi2 > CRIT_90_CHI2[len(labels)])

    excluded = np.array([x[chi2_mask] for x in raw[:-1]]).T

    store = pd.DataFrame(excluded, columns = labels)

    axes = scatter_matrix(store, alpha=0.2, figsize=(10, 10), diagonal='hist')

    if not fix_scale:
        s2 = np.log10(scale_bounds)

    if not angles:
        if fix_mixing:
            ranges = [s2, (0, 1), (0, 1), (0, 1)]
        elif fix_sfr:
            if fix_scale:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
            else:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), s2]
        else:
            if fix_scale:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
            else:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), s2, (0, 1), (0, 1), (0, 1)]
    else:
        if fix_sfr:
            if fix_mixing:
                ranges = [(0, 1), s2]
            elif fix_scale:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 2*np.pi)]
            else:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 2*np.pi), s2]
        else:
            if fix_mixing:
                ranges = [s2, (0, 1), (-1, 1)]
            elif fix_scale:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 2*np.pi), (0, 1), (-1, 1)]
            else:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 2*np.pi), s2, (0, 1), (-1, 1)]

    nof_ticks = 6
    for row_idx, row_ax in enumerate(axes):
        for col_idx, ax in enumerate(row_ax):
            xranges = ranges[col_idx]
            yranges = ranges[row_idx]
            ax.set_xlim(xranges)
            if col_idx != row_idx:
                ax.set_ylim(yranges)
            if row_idx == len(axes)-1:
                ax.tick_params(axis='x', labelsize=11)
                major_ticks = np.linspace(xranges[0], xranges[1], nof_ticks)
                ax.set_xticks(major_ticks)
                ax.xaxis.set_minor_formatter(FormatStrFormatter(""))
            if col_idx == 0:
                ax.tick_params(axis='y', labelsize=11)
                major_ticks = np.linspace(yranges[0], yranges[1], nof_ticks)[1:]
                ax.set_yticks(major_ticks)
                ax.yaxis.set_minor_formatter(FormatStrFormatter(""))

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
    # plt.figtext(0.2, 0.7, plot_label, fontsize=8)

    print 'outfile = {0}'.format(outfile)
    try:
        os.makedirs(outfile[:-len(os.path.basename(outfile))])
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(outfile[:-len(os.path.basename(outfile))]):
            pass
        else:
            raise
    plt.savefig(outfile, bbox_inches='tight', dpi=150)


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
