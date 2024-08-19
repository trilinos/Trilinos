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
# $ TRILINOS=/my/path/trilinos source stk_spack_create_env_gcc.sh

work_dir=${STK_SPACK_WORK_DIR:-/fgs/$USER/stk-spack-testing-gcc}
trilinos_source=${TRILINOS:-/fgs/$USER/Trilinos}
sierra_source=${SIERRA:-/fgs/$USER/code}

stk_spack_env=CUDA
STK_SPACK_WORK_DIR=${work_dir}
TRILINOS=${trilinos_source}
SIERRA=${sierra_source}

printf "using STK_SPACK_WORK_DIR=${STK_SPACK_WORK_DIR}\n";
printf "using TRILINOS=${TRILINOS}\n";
printf "using SIERRA=${SIERRA}\n";

if [ ! -d ${trilinos_source} ] ; then
  printf "ERROR, TRILINOS location not specified or not a directory.\n";
  return 1;
fi

if [ ! -d ${sierra_source} ] ; then
  printf "ERROR, SIERRA location not specified or not a directory.\n";
  return 1;
fi

printf "copying stk directory from SIERRA to TRILINOS...\n";
exe rm -rf ${trilinos_source}/packages/stk
exe cp -r ${sierra_source}/stk ${trilinos_source}/packages

printf "Setting up spack env 'stkSpackTesting' in STK_SPACK_WORK_DIR=${work_dir}\n"

exe mkdir -p ${work_dir}
exe cd ${work_dir}
exe rm -rf spack spack.yaml stk_test_app

exe module load aue/python/3.11.6
exe module load aue/git/2.42.0
exe module load aue/netlib-lapack/3.11.0-gcc-10.3.0
exe module load aue/openmpi/4.1.6-gcc-10.3.0

exe git clone --depth=100 --branch=releases/latest git@github.com:spack/spack.git
exe source ./spack/share/spack/setup-env.sh

exe spack env create stkSpackTesting

exe module load aue/gcc/10.3.0
exe spack compiler add
exe spack compilers

#make a copy of spack.gcc.yaml before editing it
exe cp ${sierra_source}/stk/stk_integration_tests/cmake_install_test/spack.gcc.yaml ${work_dir}/spack.yaml

spack_yaml_file=${work_dir}/spack.yaml

exe sed -i s@SED_REPLACE_INSTALL_PATH@"${work_dir}/install"@g ${spack_yaml_file}
exe sed -i s@SED_REPLACE_TRILINOS_PATH@"${trilinos_source}"@g ${spack_yaml_file}
exe spack config add -f ${spack_yaml_file}
exe spack env activate stkSpackTesting

#why do we still need the following 'spack add' commands?
#shouldn't they be loaded when we activate the environment? These specs
#are in the spack.yaml file that we just added before activating the env...

exe spack add hdf5@1.14.3~shared
exe spack add zlib
exe spack add openmpi@4.1.6
exe spack add kokkos~shared
exe spack add trilinos@develop+exodus+stk+test+gtest+kokkos~amesos~epetra~shared~boost

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

exe spack load cmake
exe spack load openmpi

printf "copying stk test app from SIERRA...\n";
exe cp -r ${sierra_source}/stk/stk_integration_tests/cmake_install_test/stk_test_app .

exe cd stk_test_app

exe source run_cmake_in_spack_env
if [ $? -ne 0 ] ; then
  printf "!! error running cmake\n";
  return 1;
fi

exe make
if [ $? -ne 0 ] ; then
  printf "!! error building\n";
  return 1;
fi

exe mpirun --np 4 ./test_stk_app
if [ $? -ne 0 ] ; then
  printf "!! error running test_stk_app\n";
  return 1;
fi

printf "all done, SUCCESS!\n";
return 0

