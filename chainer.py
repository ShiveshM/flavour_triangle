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

labels=[r'U_{e1}', r'U_{e2}', r'U_{e3}', \
        r'U_{\mu1}', r'U_{\mu2}', r'U_{\mu3}', \
        r'U_{\tau1}', r'U_{\tau2}', r'U_{\tau3}', \
        r'\phi_e', r'\phi_\mu']

ranges = []
prior_vecs = []
for x in xrange(len(labels)):
    ranges.append([0, 1])
    prior_vecs.append(np.random.uniform(0,1,10000000))

Tchain = np.load('/data/mandalia/flavour_ratio/data/mcmc_chain.npy')

Tsample=mcsamples.MCSamples(samples=Tchain,labels=labels,names=[str(i) for i in range(11)],ranges=ranges)
prior_samples=mcsamples.MCSamples(samples=prior_vecs,labels=labels,names=[str(i) for i in range(len(labels))],
                                  ranges=ranges)

g = plots.getSubplotPlotter(width_inch=10)
g.settings.legend_fontsize=20
g.settings.axes_fontsize=9
g.settings.figure_legend_frame = False
g.triangle_plot([prior_samples,Tsample], 
                filled=True, legend_labels=["Prior","Posterior"],
                legend_loc="upper right",)
g.export("fr_posterior.pdf")
