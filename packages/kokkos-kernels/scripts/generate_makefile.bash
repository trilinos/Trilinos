#!/bin/bash

KOKKOS_DEVICES=""

KOKKOS_DO_EXAMPLES="1"

KOKKOSKERNELS_OPTIONS="eti-only"
KOKKOSKERNELS_ENABLE_TPLS=""

while [[ $# > 0 ]]
do
  key="$1"

  case $key in
    --kokkoskernels-path*)
      KOKKOSKERNELS_PATH="${key#*=}"
      ;;
    --kokkos-path*)
      KOKKOS_PATH="${key#*=}"
      ;;
    --qthreads-path*)
      QTHREADS_PATH="${key#*=}"
      ;;
    --prefix*)
      PREFIX="${key#*=}"
      ;;
    --with-cuda)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Cuda"
      CUDA_PATH_NVCC=`which nvcc`
      CUDA_PATH=${CUDA_PATH_NVCC%/bin/nvcc}
      ;;
    # Catch this before '--with-cuda*'
    --with-cuda-options*)
      KOKKOS_CUDA_OPT="${key#*=}"
      ;;
    --with-cuda*)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Cuda"
      CUDA_PATH="${key#*=}"
      ;;
    --with-rocm)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},ROCm"
      ;;
    --with-openmp)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},OpenMP"
      ;;
    --with-pthread)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Pthread"
      ;;
    --with-serial)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Serial"
      ;;
    --with-qthreads*)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Qthreads"
      if [ -z "$QTHREADS_PATH" ]; then
        QTHREADS_PATH="${key#*=}"
      fi
      ;;
    --with-devices*)
      DEVICES="${key#*=}"
      KOKKOS_DEVICES="${KOKKOS_DEVICES},${DEVICES}"
      ;;
    --with-gtest*)
      GTEST_PATH="${key#*=}"
      ;;
    --with-hwloc*)
      HWLOC_PATH="${key#*=}"
      ;;
    --with-memkind*)
      KOKKOSKERNELS_SPACES="HBWSpace,${KOKKOSKERNELS_SPACES}"
      MEMKIND_PATH="${key#*=}"
      ;;
    --arch*)
      KOKKOS_ARCH="${key#*=}"
      ;;
    --cxxflags*)
      CXXFLAGS="${key#*=}"
      ;;
    --ldflags*)
      LDFLAGS="${key#*=}"
      ;;
    --debug|-dbg)
      KOKKOS_DEBUG=yes
      ;;
    --make-j*)
      echo "Warning: ${key} is deprecated"
      echo "Call make with appropriate -j flag"
      ;;
    --no-examples)
      KOKKOS_DO_EXAMPLES="0"
      ;;
    --compiler*)
      COMPILER="${key#*=}"
      CNUM=`which ${COMPILER} 2>&1 >/dev/null | grep "no ${COMPILER}" | wc -l`
      if [ ${CNUM} -gt 0 ]; then
        echo "Invalid compiler by --compiler command: '${COMPILER}'"
        exit
      fi
      if [[ ! -n  ${COMPILER} ]]; then
        echo "Empty compiler specified by --compiler command."
        exit
      fi
      CNUM=`which ${COMPILER} | grep ${COMPILER} | wc -l`
      if [ ${CNUM} -eq 0 ]; then
        echo "Invalid compiler by --compiler command: '${COMPILER}'"
        exit
      fi
      ;;
    --with-scalars*)
      KOKKOSKERNELS_SCALARS="${key#*=}"
      ;;
    --with-ordinals*)
      KOKKOSKERNELS_ORDINALS="${key#*=}"
      ;;
    --with-offsets*)
      KOKKOSKERNELS_OFFSETS="${key#*=}"
      ;;
    --with-options*)
      KOKKOSKERNELS_OPTIONS="${key#*=}"
      ;;
    --with-tpls*)
      KOKKOSKERNELS_ENABLE_TPLS="${key#*=}"
      ;;
    --with-spaces*)
      KOKKOSKERNELS_SPACES="${key#*=}"
      ;;
    --with-kokkos-options*)
      KOKKOS_OPT="${key#*=}"
      ;;
    --help)
      echo "KokkosKernels configure options:"
      echo "--kokkoskernels-path=/Path/To/KokkosKernels:  Path to the KokkosKernels root directory."
      echo "--with-scalars=[SCALARS]:             Set scalars to be instantiated."
      echo "--with-ordinals=[ORDINALS]:           Set ordinals to be instantiated."
      echo "--with-offsets=[OFFSETS]:             Set offsets to be instantiated."
      echo "--prefix=/Install/Path:               Path to install the Kokkos library."
      echo "--with-options=[OPT]:                 Set KokkosKernels Options:"
      echo "                                        eti_only: only allow ETI types to be enabled [default]"
      echo "--with-tpls=[OPT]:                    Set KokkosKernels TPLs:"
      echo "                                        mkl,blas,cublas,cusparse"
      echo ""
      echo "Kokkos configure options:"
      echo "--kokkos-path=/Path/To/Kokkos:        Path to the Kokkos root directory."
      echo "--qthreads-path=/Path/To/Qthreads:    Path to Qthreads install directory."
      echo "                                        Overrides path given by --with-qthreads."
      echo ""
      echo "--with-cuda[=/Path/To/Cuda]:          Enable Cuda and set path to Cuda Toolkit."
      echo "--with-openmp:                        Enable OpenMP backend."
      echo "--with-pthread:                       Enable Pthreads backend."
      echo "--with-serial:                        Enable Serial backend."
      echo "--with-qthreads[=/Path/To/Qthreads]:  Enable Qthreads backend."
      echo "--with-devices:                       Explicitly add a set of backends."
      echo ""
      echo "--arch=[OPT]:  Set target architectures. Options are:"
      echo "               [AMD]"
      echo "                 AMDAVX         = AMD CPU"
      echo "               [ARM]"
      echo "                 ARMv80         = ARMv8.0 Compatible CPU"
      echo "                 ARMv81         = ARMv8.1 Compatible CPU"
      echo "                 ARMv8-ThunderX = ARMv8 Cavium ThunderX CPU"
      echo "               [IBM]"
      echo "                 Power7         = IBM POWER7 and POWER7+ CPUs"
      echo "                 Power8         = IBM POWER8 CPUs"
      echo "                 Power9         = IBM POWER9 CPUs"
      echo "               [Intel]"
      echo "                 WSM            = Intel Westmere CPUs"
      echo "                 SNB            = Intel Sandy/Ivy Bridge CPUs"
      echo "                 HSW            = Intel Haswell CPUs"
      echo "                 BDW            = Intel Broadwell Xeon E-class CPUs"
      echo "                 SKX            = Intel Sky Lake Xeon E-class HPC CPUs (AVX512)"
      echo "               [Intel Xeon Phi]"
      echo "                 KNC            = Intel Knights Corner Xeon Phi"
      echo "                 KNL            = Intel Knights Landing Xeon Phi"
      echo "               [NVIDIA]"
      echo "                 Kepler30       = NVIDIA Kepler generation CC 3.0"
      echo "                 Kepler32       = NVIDIA Kepler generation CC 3.2"
      echo "                 Kepler35       = NVIDIA Kepler generation CC 3.5"
      echo "                 Kepler37       = NVIDIA Kepler generation CC 3.7"
      echo "                 Maxwell50      = NVIDIA Maxwell generation CC 5.0"
      echo "                 Maxwell52      = NVIDIA Maxwell generation CC 5.2"
      echo "                 Maxwell53      = NVIDIA Maxwell generation CC 5.3"
      echo "                 Pascal60       = NVIDIA Pascal generation CC 6.0"
      echo "                 Pascal61       = NVIDIA Pascal generation CC 6.1"
      echo ""
      echo "--compiler=/Path/To/Compiler  Set the compiler."
      echo "--debug,-dbg:                 Enable Debugging."
      echo "--cxxflags=[FLAGS]            Overwrite CXXFLAGS for library build and test"
      echo "                                build.  This will still set certain required"
      echo "                                flags via KOKKOS_CXXFLAGS (such as -fopenmp,"
      echo "                                --std=c++11, etc.)."
      echo "--ldflags=[FLAGS]             Overwrite LDFLAGS for library build and test"
      echo "                                build. This will still set certain required"
      echo "                                flags via KOKKOS_LDFLAGS (such as -fopenmp,"
      echo "                                -lpthread, etc.)."
      echo "--with-gtest=/Path/To/Gtest:  Set path to gtest.  (Used in unit and performance"
      echo "                                tests.)"
      echo "--with-hwloc=/Path/To/Hwloc:  Set path to hwloc library."
      echo "--with-memkind=/Path/To/MemKind:  Set path to memkind library."
      echo "--with-kokkos-options=[OPT]:         Additional options to Kokkos:"
      echo "                                compiler_warnings"
      echo "                                aggressive_vectorization = add ivdep on loops"
      echo "                                disable_profiling = do not compile with profiling hooks"
      echo "                                "
      echo "--with-cuda-options=[OPT]:    Additional options to CUDA:"
      echo "                                force_uvm, use_ldg, enable_lambda, rdc"
      echo "--make-j=[NUM]:               DEPRECATED: call make with appropriate"
      echo "                                -j flag"
      exit 0
      ;;
    *)
      echo "warning: ignoring unknown option $key"
      ;;
  esac

  shift
done

# Remove leading ',' from KOKKOS_DEVICES.
KOKKOS_DEVICES=$(echo $KOKKOS_DEVICES | sed 's/^,//')

# If KOKKOS_PATH undefined, assume parent dir of this script is the KOKKOS_PATH.
if [ -z "$KOKKOSKERNELS_PATH" ]; then
  KOKKOSKERNELS_PATH=$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )
else
  # Ensure KOKKOS_PATH is abs path
  KOKKOSKERNELS_PATH=$( cd $KOKKOSKERNELS_PATH && pwd )
fi

if [ "${KOKKOSKERNELS_PATH}"  = "${PWD}" ] || [ "${KOKKOSKERNELS_PATH}"  = "${PWD}/" ]; then
  echo "Running generate_makefile.sh in the KokkosKernels root directory is not allowed"
  exit
fi

KOKKOSKERNELS_SRC_PATH=${KOKKOSKERNELS_PATH}

if [ -z "$KOKKOS_PATH" ]; then
  KOKKOS_PATH=${KOKKOSKERNELS_PATH}/../kokkos
else
  # Ensure KOKKOS_PATH is abs path
  KOKKOS_PATH=$( cd $KOKKOS_PATH && pwd )
fi

#KOKKOS_SETTINGS="KOKKOS_SRC_PATH=${KOKKOS_SRC_PATH}"
#KOKKOS_SETTINGS="KOKKOS_PATH=${KOKKOS_PATH}"

if [ ${#COMPILER} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} CXX=${COMPILER}"
fi

if [ ${#KOKKOS_DEVICES} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_DEVICES=${KOKKOS_DEVICES}"
fi

if [ ${#KOKKOS_ARCH} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_ARCH=${KOKKOS_ARCH}"
fi

if [ ${#KOKKOS_DEBUG} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_DEBUG=${KOKKOS_DEBUG}"
fi

if [ ${#CUDA_PATH} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} CUDA_PATH=${CUDA_PATH}"
fi

if [ ${#CXXFLAGS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} CXXFLAGS=\"${CXXFLAGS}\""
fi

if [ ${#LDFLAGS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} LDFLAGS=\"${LDFLAGS}\""
fi

if [ ${#GTEST_PATH} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} GTEST_PATH=${GTEST_PATH}"
else
  GTEST_PATH=${KOKKOS_PATH}/tpls/gtest
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} GTEST_PATH=${GTEST_PATH}"
fi

if [ ${#HWLOC_PATH} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} HWLOC_PATH=${HWLOC_PATH}"
  KOKKOS_USE_TPLS="${KOKKOS_USE_TPLS},hwloc"
fi

if [ ${#MEMKIND_PATH} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} MEMKIND_PATH=${MEMKIND_PATH}" 
  KOKKOS_USE_TPLS="${KOKKOS_USE_TPLS},experimental_memkind"
fi

if [ ${#KOKKOS_USE_TPLS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_USE_TPLS=${KOKKOS_USE_TPLS}"
fi

if [ ${#QTHREADS_PATH} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} QTHREADS_PATH=${QTHREADS_PATH}"
fi

if [ ${#KOKKOS_OPT} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_OPTIONS=${KOKKOS_OPT}"
fi

if [ ${#KOKKOS_CUDA_OPT} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_CUDA_OPTIONS=${KOKKOS_CUDA_OPT}"
  if [[ "${KOKKOS_CUDA_OPT}" =~ "force_uvm" ]]; then
    KOKKOSKERNELS_SPACES="CudaUVMSpace,${KOKKOSKERNELS_SPACES}"
  fi
fi

if [ ${#KOKKOSKERNELS_SPACES} -gt 0 ]; then 
  KOKKOSKERNELS_SPACES="${KOKKOSKERNELS_SPACES},${KOKKOS_DEVICES}"
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOSKERNELS_SPACES=${KOKKOSKERNELS_SPACES}"  
fi

if [ ${#KOKKOSKERNELS_SCALARS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOSKERNELS_SCALARS=${KOKKOSKERNELS_SCALARS}"
fi

if [ ${#KOKKOSKERNELS_ORDINALS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOSKERNELS_ORDINALS=${KOKKOSKERNELS_ORDINALS}"
fi

if [ ${#KOKKOSKERNELS_OFFSETS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOSKERNELS_OFFSETS=${KOKKOSKERNELS_OFFSETS}"
fi

if [ ${#KOKKOSKERNELS_LAYOUTS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOSKERNELS_LAYOUTS=${KOKKOSKERNELS_LAYOUTS}"
fi

if [ ${#KOKKOSKERNELS_ENABLE_TPLS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOSKERNELS_ENABLE_TPLS=${KOKKOSKERNELS_ENABLE_TPLS}"
fi

if [ ${#KOKKOSKERNELS_OPTIONS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOSKERNELS_OPTIONS=${KOKKOSKERNELS_OPTIONS}"
fi

KOKKOS_SETTINGS_NO_KOKKOS_PATH="${KOKKOS_SETTINGS}"

KOKKOSKERNELS_TEST_INSTALL_PATH="${PWD}/install"
if [ ${#PREFIX} -gt 0 ]; then
  KOKKOSKERNELS_INSTALL_PATH="${PREFIX}"
else
  KOKKOSKERNELS_INSTALL_PATH=${KOKKOSKERNELS_TEST_INSTALL_PATH}
fi

mkdir -p install
echo "#Makefile to satisfy existens of target kokkos-clean before installing the library" > install/Makefile.kokkos
echo "kokkos-clean:" >> install/Makefile.kokkos
echo "" >> install/Makefile.kokkos
mkdir -p kokkos
mkdir -p src
mkdir -p unit_test
mkdir -p perf_test

KOKKOS_INSTALL_PATH=${PWD}/kokkos/install
KOKKOS_SETTINGS="${KOKKOS_SETTINGS_NO_KOKKOS_PATH} KOKKOS_PATH=${KOKKOS_PATH} PREFIX=${KOKKOS_INSTALL_PATH}"
echo "# KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'" > kokkos/Makefile
echo "" >> kokkos/Makefile
echo "build:" >> kokkos/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/src/Makefile ${KOKKOS_SETTINGS}" >> kokkos/Makefile
echo "" >> kokkos/Makefile
echo "install-lib:" >> kokkos/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/src/Makefile ${KOKKOS_SETTINGS} install" >> kokkos/Makefile
echo "" >> kokkos/Makefile
echo "clean:" >> kokkos/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/src/Makefile ${KOKKOS_SETTINGS} clean" >> kokkos/Makefile

KOKKOS_SETTINGS="${KOKKOS_SETTINGS_NO_KOKKOS_PATH} KOKKOS_PATH=${KOKKOS_INSTALL_PATH} KOKKOSKERNELS_PATH=${KOKKOSKERNELS_PATH} KOKKOSKERNELS_INSTALL_PATH=${KOKKOSKERNELS_INSTALL_PATH}"
# Generate subdirectory makefiles.
echo "# KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'" > src/Makefile
echo "" >> src/Makefile
echo "build:" >> src/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOSKERNELS_PATH}/src/Makefile ${KOKKOS_SETTINGS}" >> src/Makefile
echo "" >> src/Makefile
echo "install-lib:" >> src/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOSKERNELS_PATH}/src/Makefile ${KOKKOS_SETTINGS} install" >> src/Makefile
echo "" >> src/Makefile
echo "clean:" >> src/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOSKERNELS_PATH}/src/Makefile ${KOKKOS_SETTINGS} clean" >> src/Makefile

KOKKOS_SETTINGS="${KOKKOS_SETTINGS_NO_KOKKOS_PATH} KOKKOS_PATH=${KOKKOS_INSTALL_PATH} KOKKOSKERNELS_PATH=${KOKKOSKERNELS_INSTALL_PATH} KOKKOSKERNELS_SRC_PATH=${KOKKOSKERNELS_SRC_PATH}"
echo "# KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'" > unit_test/Makefile
echo "" >> unit_test/Makefile
echo "build:" >> unit_test/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOSKERNELS_PATH}/unit_test/Makefile ${KOKKOS_SETTINGS} build" >> unit_test/Makefile
echo "" >> unit_test/Makefile
echo "test: build" >> unit_test/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOSKERNELS_PATH}/unit_test/Makefile ${KOKKOS_SETTINGS} test" >> unit_test/Makefile
echo "" >> unit_test/Makefile
echo "clean:" >> unit_test/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOSKERNELS_PATH}/unit_test/Makefile ${KOKKOS_SETTINGS} clean" >> unit_test/Makefile

echo "# KOKKOS_SETTINGS='${KOKKOS_SETTINGS}'" > perf_test/Makefile
echo "" >> perf_test/Makefile
echo "build:" >> perf_test/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOSKERNELS_PATH}/perf_test/Makefile ${KOKKOS_SETTINGS} build" >> perf_test/Makefile
echo "" >> perf_test/Makefile
echo "test: build" >> perf_test/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOSKERNELS_PATH}/perf_test/Makefile ${KOKKOS_SETTINGS} test" >> perf_test/Makefile
echo "" >> perf_test/Makefile
echo "clean:" >> perf_test/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOSKERNELS_PATH}/perf_test/Makefile ${KOKKOS_SETTINGS} clean" >> perf_test/Makefile

KOKKOS_SETTINGS="${KOKKOS_SETTINGS_NO_KOKKOS_PATH} KOKKOSKERNELS_PATH=${KOKKOSKERNELS_TEST_INSTALL_PATH}"

KOKKOS_SETTINGS="${KOKKOS_SETTINGS_NO_KOKKOS_PATH} KOKKOS_PATH=${KOKKOS_PATH} KOKKOSKERNELS_PATH=${KOKKOSKERNELS_PATH}"

# Generate top level directory makefile.
echo "Generating Makefiles with options " ${KOKKOS_SETTINGS}
echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > Makefile

echo "" >> Makefile
echo "kokkos-lib:" >> Makefile
echo -e "\tcd kokkos && \$(MAKE) install-lib" >> Makefile
echo "" >> Makefile

echo "" >> Makefile
echo "kokkoskernels-lib: kokkos-lib" >> Makefile
echo -e "\tcd src && \$(MAKE) build" >> Makefile
echo "" >> Makefile

echo "install-lib: kokkoskernels-lib" >> Makefile
echo -e "\tcd src && \$(MAKE) install-lib" >> Makefile
echo "" >> Makefile

echo "build-test: install-lib" >> Makefile
echo -e "\t\$(MAKE) -C unit_test" >> Makefile
echo -e "\t\$(MAKE) -C perf_test" >> Makefile
echo "" >> Makefile

echo "test: build-test" >> Makefile
echo -e "\t\$(MAKE) -C unit_test test" >> Makefile
#echo -e "\t\$(MAKE) -C perf_test test" >> Makefile
echo "" >> Makefile

echo "clean:" >> Makefile
echo -e "\t\$(MAKE) -C unit_test clean" >> Makefile
echo -e "\t\$(MAKE) -C perf_test clean" >> Makefile
echo -e "\tcd src; \\" >> Makefile
echo -e "\t\$(MAKE) -f ${KOKKOSKERNELS_PATH}/src/Makefile ${KOKKOS_SETTINGS} clean" >> Makefile

