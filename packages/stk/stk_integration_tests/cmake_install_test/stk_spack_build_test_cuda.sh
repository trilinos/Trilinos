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
# $ SIERRA=/my/path/code source stk_spack_create_env_cuda.sh

work_dir=${STK_SPACK_WORK_DIR:-/fgs/$USER/stk-spack-testing-cuda}
sierra_source=${SIERRA:-/fgs/$USER/code}

STK_SPACK_WORK_DIR=${work_dir}
SIERRA=${sierra_source}

printf "using STK_SPACK_WORK_DIR=${STK_SPACK_WORK_DIR}\n";
printf "using SIERRA=${SIERRA}\n";

if [ ! -d ${sierra_source} ] ; then
  printf "ERROR, SIERRA location not specified or not a directory.\n";
  return 1;
fi

printf "Setting up spack env 'stkSpackTesting' in STK_SPACK_WORK_DIR=${work_dir}\n"

exe mkdir -p ${work_dir}
exe cd ${work_dir}
exe rm -rf spack spack.yaml

exe mkdir -p ${work_dir}/tmp
exe export TMPDIR=${work_dir}/tmp

exe module load aue/python/3.11.6
exe module load aue/git/2.42.0
exe module load aue/netlib-lapack/3.11.0-gcc-10.3.0
exe module load aue/openmpi/4.1.6-gcc-10.3.0

exe git clone --depth=100 git@github.com:spack/spack.git
exe source ./spack/share/spack/setup-env.sh

exe spack env create stkSpackTesting

exe module load aue/gcc/10.3.0
exe spack compiler add
exe spack compilers

#make a copy of spack.cuda.yaml before editing it
exe cp ${sierra_source}/stk/stk_integration_tests/cmake_install_test/spack.cuda.yaml ${work_dir}/spack.yaml
spack_yaml_file=${work_dir}/spack.yaml

exe sed -i s@SED_REPLACE_INSTALL_PATH@"${work_dir}/install"@g ${spack_yaml_file}
exe spack config add -f ${spack_yaml_file}
exe spack env activate stkSpackTesting

#why do we still need the following 'spack add' commands?
#shouldn't they be loaded when we activate the environment? These specs
#are in the spack.yaml file that we just added before activating the env...

exe spack add hdf5@1.14.3~shared
exe spack add zlib
exe spack add ncurses@6.3
exe spack add openmpi@4.1.6
exe spack add cuda@11.4.4
exe spack add googletest cxxstd=17
exe spack add compadre
exe spack add eigen
exe spack add yaml-cpp@0.8.0
exe spack add kokkos-kernels +cuda ~shared cuda_arch=70
exe spack add kokkos+cuda~cuda_uvm+wrapper+cuda_constexpr+cuda_lambda+cuda_relocatable_device_code~shared cuda_arch=70
exe spack add trilinos@16.1+cuda+cuda_rdc~uvm+kokkos+shards+intrepid2+stk+exodus+hdf5+zoltan2+wrapper~amesos~epetra~shared~boost cuda_arch=70 cxxstd=17

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

exe spack load googletest
exe spack load cmake
exe spack load openmpi
exe spack load yaml-cpp
exe spack load kokkos

printf "setting OMPI_CXX for CUDA environment\n";
export OMPI_CXX=$(find $(spack location -i kokkos) -name nvcc_wrapper)

printf "making build-dir for stk build...\n";
stk_build_dir=${work_dir}/build_stk
exe mkdir -p ${stk_build_dir}

exe cp ${sierra_source}/stk/stk_integration_tests/cmake_install_test/stk_test_app/run_cmake_in_spack_env ${stk_build_dir}

exe cd ${stk_build_dir}

exe STK_SOURCE_DIR=${sierra_source}/stk source run_cmake_in_spack_env

if [ $? -ne 0 ] ; then
  printf "!! error running cmake\n";
  return 1;
fi

exe make -j16
if [ $? -ne 0 ] ; then
  printf "!! error building\n";
  return 1;
fi

printf "all done, SUCCESS building!\n";
return 0

