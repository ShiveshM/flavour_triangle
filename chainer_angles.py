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

labels=[r's_{12}^2', r'c_{13}^4', r's_{23}^2', r'\delta_{CP}', r'\phi_e',
        r'\phi_\mu']
print labels

ranges = [(0, 1), (0, 1), (0, 1), (0, 2*np.pi), (0, 1), (0, 1)]


Tchain = np.load('/data/mandalia/flavour_ratio/data/mcmc_chain.npy')

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
g.export("fr_posterior_angles.pdf")
print "DONE!"
