#!/bin/bash -e

#
# This file gives the configure of all the Trilinos packages needed by Drekar
#

./attb_cuda_config_base.sh \
-D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
-D Trilinos_ENABLE_TESTS:BOOL=${TEST} \
-D Panzer_ENABLE_TESTS:BOOL=ON \
-D Trilinos_ENABLE_KokkosCore:BOOL=ON \
-D Trilinos_ENABLE_KokkosAlgorithms:BOOL=ON \
-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
-D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
-D Trilinos_ENABLE_Teko:BOOL=ON \
-D Trilinos_ENABLE_MueLu:BOOL=${EXPERIMENTAL_MUELU} \
-D Trilinos_ENABLE_Belos:BOOL=ON \
-D Trilinos_ENABLE_Panzer:BOOL=ON \
-D Trilinos_ENABLE_Shards:BOOL=ON \
-D Trilinos_ENABLE_Stratimikos:BOOL=ON \
-D Trilinos_ENABLE_ML:BOOL=ON \
-D Trilinos_ENABLE_Zoltan:BOOL=ON \
-D Trilinos_ENABLE_Zoltan2:BOOL=ON \
-D Trilinos_ENABLE_FEI:BOOL=ON \
-D Trilinos_ENABLE_Amesos:BOOL=ON \
-D Trilinos_ENABLE_SEACAS:BOOL=ON \
-D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
-D Trilinos_ENABLE_STK:BOOL=ON \
-D Trilinos_ENABLE_STKClassic:BOOL=OFF \
-D Trilinos_ENABLE_STKMesh:BOOL=ON \
-D Trilinos_ENABLE_STKUtil:BOOL=ON \
-D Trilinos_ENABLE_STKSearch:BOOL=OFF \
-D Trilinos_ENABLE_STKTopology:BOOL=ON \
-D Trilinos_ENABLE_STKTransfer:BOOL=ON \
-D Trilinos_ENABLE_STKDoc_tests:BOOL=OFF \
-D Trilinos_ENABLE_STKUnit_tests:BOOL=OFF \
-D Trilinos_ENABLE_STKUnit_test_utils:BOOL=OFF \
-D Trilinos_ENABLE_Stokhos:BOOL=OFF \
-D CMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} \
"$@"
