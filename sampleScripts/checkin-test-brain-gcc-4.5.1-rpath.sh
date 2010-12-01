#!/bin/bash

EXTRA_ARGS=$@

GCC_INSTALL_BASE=/usr/local/trilinos-toolset

echo "
-DTrilinos_EXTRA_LINK_FLAGS:STRING='-Wl,-rpath,$GCC_INSTALL_BASE/lib64' \
" > COMMON.config

echo "
-DCMAKE_CXX_COMPILER:PATH=$GCC_INSTALL_BASE/bin/g++ \ 
-DCMAKE_C_COMPILER:PATH=$GCC_INSTALL_BASE/bin/gcc \ 
" > SERIAL_RELEASE.config

../../../Trilinos/checkin-test.py \
-j12 \
--ctest-timeout=180 \
$EXTRA_ARGS  
