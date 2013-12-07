#!/bin/bash
LD_PRELOAD=/usr/lib64/libblas.so.3:\
/usr/lib64/liblapack.so.3 \
LD_LIBRARY_PATH=/home/dridzal/research/Trilinos/BUILD/gcc-4.5.1_intrepid-ML_release/install/lib:$LD_LIBRARY_PATH \
MATLABPATH=/home/dridzal/research/Trilinos/packages/intrepid/matlab/intrelab/install \
/usr/local/MATLAB/R2013b/bin/matlab -singleCompThread -nodesktop
