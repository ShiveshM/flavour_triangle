#! /usr/bin/env python
from __future__ import absolute_import, division

import numpy as np
from scipy.stats import multivariate_normal


BESTFIT=[1, 1, 1]
SIGMA=0.20

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

ls = np.linspace(0, 1, 1000)
p = cartesian((ls, ls))
print p.shape

def tc(x):
    return 1. - x[0] - x[1]
print np.array(map(tc, p))
print np.array(map(tc, p)).shape
p = np.column_stack([p, np.array(map(tc, p))])
print p.shape

p = p[p[:,-1] > 0]

cov_fr = np.identity(3) * SIGMA
llh = -np.log(multivariate_normal.pdf(p, mean=BESTFIT, cov=cov_fr))
print llh
print llh.shape

rndm = np.random.uniform(0, 65, size=len(llh))
print rndm

mask = (llh < rndm)
print mask
print np.sum(mask)

p_mask = p[mask]

with open('llh_020.txt', 'w') as f:
    for frs in p_mask:
        for fr in frs:
            f.write('{0:.6f} '.format(fr))
        f.write('\n')
