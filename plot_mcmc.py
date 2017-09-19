#! /usr/bin/env python
"""
Plot chains from MCMC
"""

from __future__ import absolute_import, division

import argparse

import numpy as np
import matplotlib
import corner

print 'Making plot'
samples = np.load('./data/mcmc_chain.npy')
fig = corner.corner(samples)
fig.savefig("untitled.png")
print 'DONE!'
