#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G

params=( "$@" )
python -u /users/mandalia/Documents/flavour_triangle/hese/mcmc_scan.py --bestfit-ratio ${params[0]} ${params[1]} ${params[2]} --sigma-ratio ${params[3]} --burnin 100 --nwalkers 200 --nsteps ${params[4]} --nufit ${params[5]} --seed 24 --outfile ${params[6]}
