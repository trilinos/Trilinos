#!/bin/bash
LD_LIBRARY_PATH=/home/dridzal/research/Trilinos/BUILD/gcc-4.5.1_intrepid-ML_release/install/lib:$LD_LIBRARY_PATH \
BLAS_VERSION=/usr/lib64/libblas.so.3 \
MATLABPATH=/home/dridzal/research/Trilinos/packages/intrepid/matlab/intrelab/install:/home/dridzal/research/Trilinos/packages/intrepid/matlab/intrelab/mesh \
/usr/local/MATLAB/R2013b/bin/matlab -singleCompThread -nodesktop
