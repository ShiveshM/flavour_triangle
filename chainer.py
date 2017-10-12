#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from getdist import plots
from getdist import mcsamples
import getdist

from matplotlib import rc, rcParams

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern'], 'size' : 18})
cols = ['#29A2C6','#FF6D31','#FFCB18','#73B66B','#EF597B', '#333333']

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 18}

labels=[r'\mid U_{e1} \mid', r'\mid U_{e2} \mid', r'\mid U_{e3} \mid', \
        r'\mid U_{\mu1} \mid', r'\mid U_{\mu2} \mid', r'\mid U_{\mu3} \mid', \
        r'\mid U_{\tau1} \mid', r'\mid U_{\tau2} \mid', r'\mid U_{\tau3} \mid', \
        r'\phi_e', r'\phi_\mu']
print labels

ranges = []
for x in xrange(len(labels)):
    ranges.append([0, 1])

Tchain = np.load('/data/mandalia/flavour_ratio/data/mcmc_chain.npy')
# Tchain = np.load('./data/mcmc_chain.npy')

Tsample=mcsamples.MCSamples(
    samples=Tchain, labels=labels, ranges=ranges
)

Tsample.updateSettings({'contours': [0.90, 0.99]})
Tsample.num_bins_2D=500
Tsample.fine_bins_2D=500
Tsample.smooth_scale_2D=0.03

g = plots.getSubplotPlotter()
g.settings.num_plot_contours = 2
g.settings.legend_fontsize = 25
g.settings.axes_fontsize = 10
g.settings.figure_legend_frame = False
g.triangle_plot(
    [Tsample], filled=True, legend_loc="upper right"
)
g.export("fr_posterior.pdf")
