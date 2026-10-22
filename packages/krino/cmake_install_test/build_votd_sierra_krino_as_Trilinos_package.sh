#!/bin/sh

sierra_proj=${SIERRA_PROJ:-${PWD}}
output_dir=${OUTPUT_DIR:-${PWD}/../krino-cmake-testing}

source $sierra_proj/krino/cmake_install_test/build_helpers.sh

cuda_on_or_off=${CUDA:-OFF}
date_suffix=`date +%F_%H-%M-%S`

setup_build_with_krino
setup_environment

build_yaml release

enable_yaml=ON
enable_pykrino=ON
build_trilinos_with_krino release
run_tests release
build_trilinos_with_krino debug
run_tests debug

enable_yaml=OFF
enable_pykrino=OFF
build_trilinos_with_krino release
run_tests release
