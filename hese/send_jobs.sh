#!/bin/bash
unset module;

a_ifr=(1 2 0)
b_ifr=(1 0 0)
c_ifr=(0 1 0)
d_ifr=(0 0 1)
m_a=(
a_ifr[@]
b_ifr[@]
c_ifr[@]
d_ifr[@]
)

sigma=0.01
nsteps=2000
outfile=/data/mandalia/flavour_ratio/data/mcmc_chain

count=${#m_a[@]}
for ((i=0; i<$count; i++)); do
    arg=${!m_a[i]}
    echo "${arg[0]} ${arg[1]} ${arg[2]}"

    # /users/mandalia/Documents/flavour_triangle/hese/wrap.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} ${nsteps} False ${outfile}
    # break

    qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/wrap.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} ${nsteps} False ${outfile}
    qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/wrap.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} ${nsteps} True ${outfile}
    # qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/wrap_ifr.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} ${arg[0]} ${arg[1]} ${arg[2]} ${nsteps} False ${outfile}
    # qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/wrap_ifr.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} ${arg[0]} ${arg[1]} ${arg[2]} ${nsteps} True ${outfile}

    # for brg in a_ifr b_ifr c_ifr d_ifr; do
	# /users/mandalia/Documents/flavour_triangle/hese/wrap_ifr.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} ${brg[0]} ${brg[1]} ${brg[2]} ${nsteps} False ${outfile}
	# /users/mandalia/Documents/flavour_triangle/hese/wrap_ifr.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} ${brg[0]} ${brg[1]} ${brg[2]} ${nsteps} True ${outfile}
    # done

    # break
done
