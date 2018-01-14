#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=30G

params=( "$@" )
python -u /users/mandalia/Documents/flavour_triangle/hese/bsm/mcmc_scan_afr.py --bestfit-ratio ${params[0]} ${params[1]} ${params[2]} --sigma-ratio ${params[3]} --burnin 20000 --nwalkers 200 --nsteps ${params[4]} --nufit ${params[5]} --seed 24 --outfile ${params[6]}
