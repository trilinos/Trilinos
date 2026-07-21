#!/usr/bin/bash

exe() {
  stdbuf -o0 -e0 echo "% $@" ;
  eval "$@" ;
  if [ $? -ne 0 ] ; then
    echo "'$@' failed.";
    return 1;
  fi
}

# To specify custom paths for one or more of the following, set
# the variable on the command line when running this script.
# Example:
# $ TRILINOS=/my/path/trilinos source stk_spack_create_env.sh

work_dir=${STK_SPACK_WORK_DIR:-/fgs/$USER/stk-cmake-testing}
sierra_dir=${SIERRA_DIR:-/fgs/$USER/code}
skip_spack_setup=${SKIP_SPACK_SETUP:-false}

STK_SPACK_WORK_DIR=${work_dir}

printf "using STK_SPACK_WORK_DIR=${STK_SPACK_WORK_DIR}\n";

exe cd ${STK_SPACK_WORK_DIR}

if [ "${skip_spack_setup}" == "false" ] ; then
exe rm -rf spack
exe git clone --depth=100 git@github.com:spack/spack.git
fi

exe source ./spack/share/spack/setup-env.sh

if [ "${skip_spack_setup}" == "false" ] ; then
  exe spack env create stk_gcc12_spack_env
fi

exe spack env activate stk_gcc12_spack_env

exe cp ${sierra_dir}/stk/stk_integration_tests/cmake_install_test/spack.gcc.yaml spack.yaml
exe sed -i s@SED_REPLACE_STK_SPACK_INSTALL_PATH@"${work_dir}/spack-install"@g spack.yaml
exe sed -i s@SED_REPLACE_STK_SPACK_STAGE_PATH@"${work_dir}/spack-stage"@g spack.yaml
exe spack config add -f spack.yaml

exe spack add gcc@12.3.0
exe spack add cmake@3.29.4
exe spack add boost@1.85.0
exe spack add googletest cxxstd=20
exe spack add zlib@1.3
exe spack add openmpi@4.1.6
exe spack add hdf5@1.14.3
exe spack add parallel-netcdf@1.12.3
exe spack add netcdf-c@4.9.2
exe spack add metis@5.1.0
exe spack add parmetis@4.0.3
exe spack add yaml-cpp@0.8.0

exe spack concretize -f

if [ $? -ne 0 ] ; then
  printf "!! error running spack concretize\n";
  return 1;
fi

exe spack install

if [ $? -ne 0 ] ; then
  printf "!! error running spack install\n";
  return 1;
fi

exe spack load gcc
exe spack load cmake
exe spack load openmpi


printf "all done, SUCCESS!\n";
return 0

