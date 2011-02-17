#!/bin/bash

# Used to test Trilinos on u233.ornl.gov

EXTRA_ARGS=$@

echo "
-DBUILD_SHARED_LIBS:BOOL=ON
-DTPL_BLAS_LIBRARIES=/usr/lib64/libblas.so.3
-DTPL_LAPACK_LIBRARIES=/usr/lib64/liblapack.so.3
" > COMMON.config

../../Trilinos/checkin-test.py \
-j12 \
--ctest-timeout=180 \
--ctest-options="-E Teuchos_Dependencies_Serialization_test_MPI_1" \
$EXTRA_ARGS  
