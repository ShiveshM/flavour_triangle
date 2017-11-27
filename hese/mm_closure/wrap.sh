#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G

params=( "$@" )
python -u /users/mandalia/Documents/flavour_triangle/hese/mm_closure/mcmc_scan.py --injected-ratio ${params[0]} ${params[1]} --sigma-ratio ${params[2]} --hypothesis-ratio ${params[3]} ${params[4]} --burnin 1000 --nwalkers 200 --nsteps ${params[5]} --nufit False --seed 24 --outfile ${params[6]}
