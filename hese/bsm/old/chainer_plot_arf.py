#! /usr/bin/env python
"""
From an MCMC chains file, make a triangle plot.
"""

from __future__ import absolute_import, division

import sys
sys.path.append('/users/mandalia/Documents/flavour_triangle/hese/bsm/')

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import rc, rcParams

import getdist
from getdist import plots
from getdist import mcsamples

import mcmc_scan_afr

DIMENSION = 3
ENERGY = 1000 # GeV
SCALE_BOUNDS = [(1e-25, 1e-19), (1e-23, 1e-17)]

rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})
cols = ['#29A2C6','#FF6D31','#FFCB18','#73B66B','#EF597B', '#333333']

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 18}


def plot(infile, angles, nufit, outfile, bestfit_ratio=None, sigma_ratio=None):
    """Make the triangle plot"""
    if not angles:
        labels = [r'\mid \tilde{U}_{e1} \mid', r'\mid \tilde{U}_{e2} \mid', r'\mid \tilde{U}_{e3} \mid', \
                  r'\mid \tilde{U}_{\mu1} \mid', r'\mid \tilde{U}_{\mu2} \mid', r'\mid \tilde{U}_{\mu3} \mid', \
                  r'\mid \tilde{U}_{\tau1} \mid', r'\mid \tilde{U}_{\tau2} \mid', r'\mid \tilde{U}_{\tau3} \mid', \
                  r'\Lambda_1', r'\Lambda_2', r'\phi_e', r'\phi_\mu', r'\phi_\tau']
    else:
        labels=[r's_{12}^2', r'c_{13}^4', r's_{23}^2', r'\delta_{CP}',
                r'\Lambda_1', r'\Lambda_2', r'sin^4(\phi)', r'cos(2\psi)']
    print labels

    s1 = np.log10(SCALE_BOUNDS[0])
    s2 = np.log10(SCALE_BOUNDS[1])
    if not angles and nufit:
        assert 0
        ranges = [(0.800, 0.844), (0.515, 0.581), (0.139, 0.155), (0.229, 0.516),
                  (0.438, 0.699), (0.614, 0.790), (0.249, 0.528), (0.462, 0.715),
                  # (0.595, 0.776), (0, 1), (0, 1)]
                  (0.595, 0.776), (0, 1), (0, 1), (0, 1)]
    elif not angles:
        ranges = [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), s1, s2, (0, 1), (0, 1), (0, 1)]
    else:
        ranges = [(0, 1), (0, 1), (0, 1), (0, 2*np.pi), s1, s2, (0, 1), (-1, 1)]

    def flat_angles_to_u(x):
        return abs(mcmc_scan_afr.angles_to_u(x)).astype(np.float32).flatten().tolist()

    raw = np.load(infile)
    if not angles:
        fr_elements = np.array(map(mcmc_scan_afr.angles_to_fr, raw[:,-2:]))
	sc_elements = raw[:,-4:-2]
        m_elements = np.array(map(flat_angles_to_u, raw[:,:-4]))
        Tchain = np.column_stack([m_elements, sc_elements, fr_elements])
    else:
        Tchain = raw

    if bestfit_ratio is not None and sigma_ratio is not None:
        label = 'Bestfit ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nSigma = {3:.3f}\nDimension = {4}\nEnergy = {5} GeV'.format(
            bestfit_ratio[0], bestfit_ratio[1], bestfit_ratio[2], sigma_ratio, DIMENSION, ENERGY
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
        '--nufit', default=False, action='store_true',
        help='Include NuFit priors'
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
