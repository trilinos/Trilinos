#!/bin/bash

# ./testboth.sh > output.txt
# grep ">>>\|>> \|Timer:" output.txt

#numacmd="KMP_AFFINITY=balanced numactl --membind 1"

# 128x128x128 domain is too large to use hbm only
numacmd="KMP_AFFINITY=balanced"
#sz="-ni 128 -nj 128 -nk 128"

# 4 5 8 9 A < 10 GB
for bsz in 3 5; do
    sz="-ni 128 -nj 128 -nk 128"
    echo ">>> bsz $bsz"
    for nth in 4 8 16 34 68 136 272; do
        echo ">> nthread $nth"
        echo "> kk"
        cmd="$numacmd ./KokkosKernels_Test_BlockCrs --kokkos-threads=$nth $sz -bs $bsz"
        echo $cmd
        eval $cmd
        echo "> sparc"
        cmd="$numacmd ./bcrs --kokkos-threads=$nth $sz -bs $bsz"
        echo $cmd
        eval $cmd
    done
done

for bsz in 10 15; do
    sz="-ni 64 -nj 64 -nk 128"        
    echo ">>> bsz $bsz"
    for nth in 4 8 16 34 68 136 272; do
        echo ">> nthread $nth"
        echo "> kk"
        cmd="$numacmd ./KokkosKernels_Test_BlockCrs --kokkos-threads=$nth $sz -bs $bsz"
        echo $cmd
        eval $cmd
        echo "> sparc"
        cmd="$numacmd ./bcrs --kokkos-threads=$nth $sz -bs $bsz"
        echo $cmd
        eval $cmd
    done
done
