#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G

params=( "$@" )
python -u /users/mandalia/Documents/flavour_triangle/hese/mcmc_scan_ifr.py --bestfit-ratio ${params[0]} ${params[1]} ${params[2]} --sigma-ratio ${params[3]} --hypothesis-ratio ${params[4]} ${params[5]} ${params[6]} --burnin 4000 --nwalkers 200 --nsteps ${params[7]} --nufit ${params[8]} --seed 24 --outfile ${params[9]}
