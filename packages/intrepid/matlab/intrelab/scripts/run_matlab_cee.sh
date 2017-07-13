#!/bin/bash
LD_PRELOAD=/projects/sems/install/rhel6-x86_64/sems/compiler/gcc/5.3.0/base/lib64/libstdc++.so: \
LD_LIBRARY_PATH=/scratch/dridzal/Trilinos/BUILD/intrelab/install/lib:$LD_LIBRARY_PATH \
BLAS_VERSION=/usr/lib64/libblas.so \
MATLABPATH=/scratch/dridzal/Trilinos/packages/intrepid/matlab/intrelab/install:/scratch/dridzal/Trilinos/packages/intrepid/matlab/intrelab/mesh \
/usr/local/matlab/8.6b/bin/matlab -singleCompThread -nodesktop -nosplash
