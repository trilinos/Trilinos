#!/bin/bash
LD_PRELOAD=/projects/sems/install/rhel6-x86_64/sems/compiler/gcc/5.3.0/base/lib64/libstdc++.so: \
LD_LIBRARY_PATH=/home/dridzal/research/Trilinos/BUILD/intrelab/install/lib:$LD_LIBRARY_PATH \
BLAS_VERSION=/usr/lib64/libblas.so \
MATLABPATH=/home/dridzal/research/Trilinos/packages/intrepid/matlab/intrelab/install:/home/dridzal/research/Trilinos/packages/intrepid/matlab/intrelab/mesh \
/home/dridzal/software/MATLAB/R2016b/bin/matlab -singleCompThread -nodesktop -nosplash

