#! /usr/bin/env python
from __future__ import absolute_import, division

import numpy as np
from scipy.stats import multivariate_normal


BESTFIT=[1, 1, 1]
SIGMA=0.2

def cartesian(arrays, out=None):
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

ls = np.linspace(0, 1, 100)
p = cartesian((ls, ls, ls))
print p.shape

cov_fr = np.identity(3) * SIGMA
llh = -np.log10(multivariate_normal.pdf(p, mean=BESTFIT, cov=cov_fr))
print llh
print llh.shape

rndm = np.random.uniform(np.min(llh), np.max(llh), size=len(llh))
print rndm

mask = (llh > rndm)
print mask
print np.sum(mask)

p_mask = p[mask]

with open('llh.txt', 'w') as f:
    for frs in p_mask:
        for fr in frs:
            f.write('{0:.6f} '.format(fr))
        f.write('\n')
