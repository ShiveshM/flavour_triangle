#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=30G

params=( "$@" )
python -u /users/mandalia/Documents/flavour_triangle/hese/bsm/mcmc_scan_ifr.py --bestfit-ratio ${params[0]} ${params[1]} ${params[2]} --sigma-ratio ${params[3]} --hypothesis-ratio ${params[4]} ${params[5]} ${params[6]} --burnin 2000 --nwalkers 200 --nsteps ${params[7]} --nufit ${params[8]} --seed 24 --outfile ${params[9]}
