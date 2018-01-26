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
# g_ifr[@]
e_ifr[@]
# f_ifr[@]
# a_ifr[@]
# b_ifr[@]
c_ifr[@]
# d_ifr[@]
)

com_a=(1 1 1 1 0 0)
com_b=(1 1 1 0 1 0)
com_c=(1 1 1 0 0 1)
com_d=(1 1 1 1 2 0)
com_e=(1 1 0 0 1 0)
com_f=(1 1 0 1 2 0)
com_g=(1 1 0 1 0 0)
com_h=(1 0 0 1 0 0)
com_i=(0 1 0 0 1 0)
com_j=(1 2 0 0 1 0)
com_k=(1 2 0 1 2 0)

a_c=(
com_a[@]
com_b[@]
com_c[@]
com_d[@]
com_e[@]
com_f[@]
com_g[@]
com_h[@]
com_i[@]
com_j[@]
com_k[@]
)

sigma=0.1
dimension=3
energy=100000
flat=False
burnin=2000
nwalkers=200
nsteps=200000
#scales="1E-20 1E-23 1E-25 1E-27 1E-30"
# scales="5E-24 1E-24 5E-25 1E-25 5E-26 1E-26 5E-27 1E-27 5E-28 1E-28 5E-29 1E-29 5E-30 1E-30"
scales="1E-20 1E-30"

count=${#a_c[@]}
for ((i=0; i<$count; i++));do
    arg=(${!a_c[i]})
    echo "${arg[0]} ${arg[1]} ${arg[2]} ${arg[3]} ${arg[4]} ${arg[5]}"
    outfile=/data/mandalia/flavour_ratio/data/DIM${dimension}/fix_ifr/0_1/mcmc_chain
    qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/bsm/wrap.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} True ${arg[3]} ${arg[4]} ${arg[5]} False 0 ${dimension} ${energy} ${flat} ${burnin} ${nwalkers} ${nsteps} ${outfile} False
done

count=${#m_a[@]}
for ((i=0; i<$count; i++)); do
    arg=(${!m_a[i]})
    echo "${arg[0]} ${arg[1]} ${arg[2]}"

    outfile=/data/mandalia/flavour_ratio/data/DIM${dimension}/full_scan/0_1/mcmc_chain
    qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/bsm/wrap.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} False 0 0 0 False 0 ${dimension} ${energy} ${flat} ${burnin} ${nwalkers} ${nsteps} ${outfile} False

    # outfile=/data/mandalia/flavour_ratio/data/DIM${dimension}/fix_mixing_zeromass/0_001/mcmc_chain
    # qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/bsm/wrap.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} False 0 0 0 False 0 ${dimension} ${energy} ${flat} ${burnin} ${nwalkers} ${nsteps} ${outfile} True

    # for ((j=0; j<$count; j++)); do
    #     brg=${!m_a[j]}
# #         # brg=(${!m_a[i]})
    #     echo "arg = $arg , brg = $brg"

        # outfile=/data/mandalia/flavour_ratio/data/DIM${dimension}/fix_ifr/0_001/mcmc_chain
        # qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/bsm/wrap.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} True ${brg[0]} ${brg[1]} ${brg[2]} False 0 ${dimension} ${energy} ${flat} ${burnin} ${nwalkers} ${nsteps} ${outfile} False

#         # outfile=/data/mandalia/flavour_ratio/data/DIM${dimension}/fix_ifr_mixing/0_001/mcmc_chain
#         # qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/bsm/wrap.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} True ${brg[0]} ${brg[1]} ${brg[2]} False 0 ${dimension} ${energy} ${flat} ${burnin} ${nwalkers} ${nsteps} ${outfile} True

# 	# for scale in ${scales}; do
#             # echo "${scale}"
#             # outfile=/data/mandalia/flavour_ratio/data/DIM${dimension}/fix_ifr_scale/0_001/mcmc_chain
# 	#     qsub -cwd -V /users/mandalia/Documents/flavour_triangle/hese/bsm/wrap.sh ${arg[0]} ${arg[1]} ${arg[2]} ${sigma} True ${brg[0]} ${brg[1]} ${brg[2]} True ${scale} ${dimension} ${energy} ${flat} ${burnin} ${nwalkers} ${nsteps} ${outfile} False
# 	# done
        # break
    # done

    # break
done
