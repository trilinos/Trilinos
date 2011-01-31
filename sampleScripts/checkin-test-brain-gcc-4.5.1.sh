#!/bin/bash

EXTRA_ARGS=$@

echo "
" > COMMON.config

echo "
-DCMAKE_C_COMPILER:PATH=/usr/local/gcc-4.5.1/bin/gcc
-DCMAKE_CXX_COMPILER:PATH=/usr/local/gcc-4.5.1/bin/g++
" > SERIAL_RELEASE.config

export LD_LIBRARY_PATH=/usr/local/gcc-4.5.1/lib64 

../../../Trilinos/checkin-test.py \
-j12 \
--ctest-timeout=180 \
$EXTRA_ARGS  
