#!/bin/bash
LD_PRELOAD=/home/dridzal/software/gcc/install_gcc-4.5.1/lib64/libstdc++.so.6:\
/home/dridzal/software/gcc/install_gcc-4.5.1/lib64/libgcc_s.so.1:\
/home/dridzal/software/gcc/install_gcc-4.5.1/lib64/libgfortran.so.3:\
/usr/lib64/libblas.so.3:\
/usr/lib64/liblapack.so.3 \
LD_LIBRARY_PATH=/home/dridzal/research/Trilinos/BUILD/gcc-4.5.1_intrepid-ML_release/install/lib:$LD_LIBRARY_PATH \
MATLABPATH=/home/dridzal/research/Trilinos/packages/intrepid/matlab/intrelab/install:/home/dridzal/research/Trilinos/packages/intrepid/matlab/intrelab/mesh \
/usr/local/MATLAB/R2010b/bin/matlab -singleCompThread -nodesktop
