#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G

params=( "$@" )
python -u /users/mandalia/Documents/flavour_triangle/hese/test_haar/mcmc_scan.py --bestfit-ratio ${params[0]} ${params[1]} ${params[2]} --sigma-ratio ${params[3]} --burnin 1000 --nwalkers 200 --nsteps ${params[4]} --seed 24 --outfile ${params[5]}
