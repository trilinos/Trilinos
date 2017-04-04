#!/bin/sh -x

module purge
module load sierra-devel

cp $1 bridge.png

/scratch/sierra/code/bin/linux-gcc-4.9.3-ip-openmpi-1.6.4/release/stk_learning_game_of_life_utest --gtest_filter=TOSDTWD.hex_mesh_from_image_multiple_blocks -i bridge.png

module purge
module load sierra

sierra -j 6 adagio -i bridgeIt.i

module load viz
ens -p lookAtMe.cmd.enc
