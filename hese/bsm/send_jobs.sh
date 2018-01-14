#!/bin/bash
unset module;

a_ifr=(1 2 0)
b_ifr=(1 0 0)
c_ifr=(0 1 0)
d_ifr=(0 0 1)
e_ifr=(1 1 1)
f_ifr=(2 1 0)
g_ifr=(1 1 0)

# m_a=(
# a_ifr[@]
# b_ifr[@]
# c_ifr[@]
# d_ifr[@]
# )
m_a=(
g_ifr[@]
# e_ifr[@]
# f_ifr[@]
# a_ifr[@]
# b_ifr[@]
c_ifr[@]
# d_ifr[@]
)

sigma=0.001
dimension=3
energy=100000
flat=False
burnin=200
nwalkers=200
nsteps=2000
scales="1E-20 1E-25 1E-30"

count=${#m_a[@]}
for ((i=0; i<$count; i++)); do
    arg=${!m_a[i]}
    echo "${arg[0]} ${arg[1]} ${arg[2]}"

    outfile=/data/mandalia/flavour_ratio/data/DIM${dimension}/full_scan/0_001/mcmc_chain
    qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/bsm/wrap.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} False 0 0 0 False 0 ${dimension} ${energy} ${flat} ${burnin} ${nwalkers} ${nsteps} ${outfile}

    for ((j=0; j<$count; j++)); do
        brg=${!m_a[j]}

        outfile=/data/mandalia/flavour_ratio/data/DIM${dimension}/fix_ifr/0_001/mcmc_chain
        qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/bsm/wrap.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} True ${brg[0]} ${brg[1]} ${brg[2]} False 0 ${dimension} ${energy} ${flat} ${burnin} ${nwalkers} ${nsteps} ${outfile}

	for scale in ${scales}; do
            echo "${scale}"
            outfile=/data/mandalia/flavour_ratio/data/DIM${dimension}/fix_ifr_scale/0_001/mcmc_chain
	    qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/bsm/wrap.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} True ${brg[0]} ${brg[1]} ${brg[2]} True ${scale} ${dimension} ${energy} ${flat} ${burnin} ${nwalkers} ${nsteps} ${outfile}
	done

        # OLD
	# qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/bsm/wrap_ifr.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} ${brg[0]} ${brg[1]} ${brg[2]} ${nsteps} False ${outfile}
	# /users/mandalia/Documents/flavour_triangle/hese/bsm/wrap_ifr.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} ${brg[0]} ${brg[1]} ${brg[2]} ${nsteps} False ${outfile}
    done

    # OLD
    # qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/bsm/wrap_afr.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} ${nsteps} False ${outfile}

    # break
done
