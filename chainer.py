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
# labels=[r's_{12}^2', r'c_{13}^4', r's_{23}^2', r'\delta_{CP}', r'\phi_e',
#         r'\phi_\mu']
print labels

# ranges = [(0, 1), (0, 1), (0, 1), (0, 2*np.pi), (0, 1), (0, 1)]
ranges = []
for x in xrange(len(labels)):
    ranges.append([0, 1])


def angles_to_u(angles):
    s12_2, c13_4, s23_2, dcp = angles
    dcp = np.complex128(dcp)

    c13_2 = np.sqrt(c13_4)

    c12 = np.sqrt(1. - s12_2)
    s12 = np.sqrt(s12_2)
    c13 = np.sqrt(c13_2)
    s13 = np.sqrt(1. - c13_2)
    c23 = np.sqrt(1. - s23_2)
    s23 = np.sqrt(s23_2)

    phase = np.exp(-1*dcp)
    p1 = np.array([[1   , 0   , 0]         , [0    , c23 , s23] , [0          , -s23 , c23]])
    p2 = np.array([[c13 , 0   , s13*phase] , [0    , 1   , 0]   , [-s13*phase , 0    , c13]])
    p3 = np.array([[c12 , s12 , 0]         , [-s12 , c12 , 0]   , [0          , 0    , 1]])

    # u = np.dot(p1, p2, p3)
    u = abs(np.dot(np.dot(p1, p2), p3)).astype(np.float32)
    return u.flatten().tolist()


raw = np.load('/data/mandalia/flavour_ratio/data/mcmc_chain.npy')
angles = np.array(map(angles_to_u, raw[:,:-2]))
Tchain = np.hstack([angles, raw[:,-2:]])

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
print "DONE!"
