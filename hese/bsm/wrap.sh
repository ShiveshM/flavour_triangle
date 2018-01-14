#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=30G

params=( "$@" )
echo "python -u /users/mandalia/Documents/flavour_triangle/hese/bsm/mcmc_scan.py --measured-ratio ${params[0]} ${params[1]} ${params[2]} --sigma-ratio ${params[3]} --fix-source-ratio ${params[4]} --source-ratio ${params[5]} ${params[6]} ${params[7]} --fix-scale ${params[8]} --scale ${params[9]} --dimension ${params[10]} --energy ${params[11]} --flat-llh ${params[12]} --burnin ${params[13]} --nwalkers ${params[14]} --nsteps ${params[15]} --seed 24 --outfile ${params[16]}"
python -u /users/mandalia/Documents/flavour_triangle/hese/bsm/mcmc_scan.py --measured-ratio ${params[0]} ${params[1]} ${params[2]} --sigma-ratio ${params[3]} --fix-source-ratio ${params[4]} --source-ratio ${params[5]} ${params[6]} ${params[7]} --fix-scale ${params[8]} --scale ${params[9]} --dimension ${params[10]} --energy ${params[11]} --flat-llh ${params[12]} --burnin ${params[13]} --nwalkers ${params[14]} --nsteps ${params[15]} --seed 24 --outfile ${params[16]}
