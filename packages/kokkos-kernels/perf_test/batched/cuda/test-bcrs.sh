#!/bin/bash

# ./testboth.sh > output.txt
# grep ">>>\|>> \|Timer:" output.txt

# 128x128x128 domain is too large to use hbm only
#sz="-ni 128 -nj 128 -nk 128"

# 4 5 8 9 A < 10 GB
for bsz in 3 5; do
    sz="-ni 128 -nj 128 -nk 128"
    echo ">>> bsz $bsz"

    echo "> kk"
    cmd="./KokkosKernels_Test_BlockCrs $sz -bs $bsz -opf 1 -ops 1"
    echo $cmd
    eval $cmd 

    echo "> sparc"
    cmd="./bcrs $sz -bs $bsz"
    echo $cmd
    eval $cmd
done


for bsz in 10 15; do
    sz="-ni 64 -nj 64 -nk 128"        
    echo ">>> bsz $bsz"

    echo "> kk"
    cmd="./KokkosKernels_Test_BlockCrs $sz -bs $bsz -opf 1 -ops 1"
    echo $cmd
    eval $cmd

    echo "> sparc"
    cmd="./bcrs $sz -bs $bsz"
    echo $cmd
    eval $cmd
done
