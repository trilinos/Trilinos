#!/bin/bash

update_kokkos_devices() {
   SEARCH_TEXT="*$1*"
   if [[ $KOKKOS_DEVICES == $SEARCH_TEXT ]]; then
      echo kokkos devices already includes $SEARCH_TEXT
   else
      if [ "$KOKKOS_DEVICES" = "" ]; then
         KOKKOS_DEVICES="$1"
         echo reseting kokkos devices to $KOKKOS_DEVICES
      else
         KOKKOS_DEVICES="${KOKKOS_DEVICES},$1"
         echo appending to kokkos devices $KOKKOS_DEVICES
      fi
   fi
}

get_kokkos_device_list() {
  KOKKOS_DEVICE_CMD=
  PARSE_DEVICES_LST=$(echo $KOKKOS_DEVICES | tr "," "\n")
  for DEVICE_ in $PARSE_DEVICES_LST
  do 
     UC_DEVICE=$(echo $DEVICE_ | tr "[:lower:]" "[:upper:]")
     KOKKOS_DEVICE_CMD="-DKokkos_ENABLE_${UC_DEVICE}=ON ${KOKKOS_DEVICE_CMD}"
  done
}

get_kokkos_arch_list() {
  KOKKOS_ARCH_CMD=
  PARSE_ARCH_LST=$(echo $KOKKOS_ARCH | tr "," "\n")
  for ARCH_ in $PARSE_ARCH_LST
  do 
     UC_ARCH=$(echo $ARCH_ | tr "[:lower:]" "[:upper:]")
     KOKKOS_ARCH_CMD="-DKokkos_ARCH_${UC_ARCH}=ON ${KOKKOS_ARCH_CMD}"
  done
}

get_kokkos_cuda_option_list() {
  echo parsing KOKKOS_CUDA_OPTIONS=$KOKKOS_CUDA_OPTIONS
  KOKKOS_CUDA_OPTION_CMD=
  PARSE_CUDA_LST=$(echo $KOKKOS_CUDA_OPTIONS | tr "," "\n")
  for CUDA_ in $PARSE_CUDA_LST
  do 
     CUDA_OPT_NAME=
     if [ "${CUDA_}" == "enable_lambda" ]; then
        CUDA_OPT_NAME=CUDA_LAMBDA
     elif  [ "${CUDA_}" == "rdc" ]; then	
        CUDA_OPT_NAME=CUDA_RELOCATABLE_DEVICE_CODE
     elif  [ "${CUDA_}" == "force_uvm" ]; then
        CUDA_OPT_NAME=CUDA_UVM
     elif  [ "${CUDA_}" == "use_ldg" ]; then
        CUDA_OPT_NAME=CUDA_LDG_INTRINSIC
     else
        echo "${CUDA_} is not a valid cuda options..."
     fi
     if [ "${CUDA_OPT_NAME}" != "" ]; then
        KOKKOS_CUDA_OPTION_CMD="-DKokkos_ENABLE_${CUDA_OPT_NAME}=ON ${KOKKOS_CUDA_OPTION_CMD}"
     fi
  done
}

get_kokkos_option_list() {
  echo parsing KOKKOS_OPTIONS=$KOKKOS_OPTIONS
  KOKKOS_OPTION_CMD=
  PARSE_OPTIONS_LST=$(echo $KOKKOS_OPTIONS | tr "," "\n")
  for OPT_ in $PARSE_OPTIONS_LST
  do 
     UC_OPT_=$(echo $OPT_ | tr "[:lower:]" "[:upper:]")
     if [[ "$UC_OPT_" == *DISABLE* ]]; then
        FLIP_OPT_=${UC_OPT_/DISABLE/ENABLE}
        KOKKOS_OPTION_CMD="-DKokkos_${FLIP_OPT_}=OFF ${KOKKOS_OPTION_CMD}"
     elif [[ "$UC_OPT_" == *ENABLE* ]]; then
        KOKKOS_OPTION_CMD="-DKokkos_${UC_OPT_}=ON ${KOKKOS_OPTION_CMD}"
     else
        KOKKOS_OPTION_CMD="-DKokkos_ENABLE_${UC_OPT_}=ON ${KOKKOS_OPTION_CMD}"
     fi
  done
}

get_kernels_scalar_list() {
  echo "parsing KOKKOSKERNELS_SCALARS=$KOKKOSKERNELS_SCALARS"
  KOKKOSKERNELS_SCALARS_CMD=
  PARSE_SCALARS_LIST=$(echo $KOKKOSKERNELS_SCALARS | tr "," "\n")
  for SCALAR_ in $PARSE_SCALARS_LIST
  do
    UC_SCALAR=$(echo $SCALAR_ | tr "[:lower:]" "[:upper:]")
    KOKKOSKERNELS_SCALARS_CMD="-DKokkosKernels_INST_${UC_SCALAR}=ON ${KOKKOSKERNELS_SCALARS_CMD}"
  done
}

get_kernels_ordinals_list() {
  echo "parsing KOKKOSKERNELS_ORDINALS=$KOKKOSKERNELS_ORDINALS"
  KOKKOSKERNELS_ORDINALS_CMD=
  PARSE_ORDINALS_LIST=$(echo $KOKKOSKERNELS_ORDINALS | tr "," "\n")
  for ORDINALS_ in $PARSE_ORDINALS_LIST
  do
    UC_ORDINALS=$(echo $ORDINALS_ | tr "[:lower:]" "[:upper:]")
    KOKKOSKERNELS_ORDINALS_CMD="-DKokkosKernels_INST_ORDINAL_${UC_ORDINALS}=ON ${KOKKOSKERNELS_ORDINALS_CMD}"
  done
}

get_kernels_offsets_list() {
  echo "parsing KOKKOSKERNELS_OFFSETS=$KOKKOSKERNELS_OFFSETS"
  KOKKOSKERNELS_OFFSETS_CMD=
  PARSE_OFFSETS_LIST=$(echo $KOKKOSKERNELS_OFFSETS | tr "," "\n")
  for OFFSETS_ in $PARSE_OFFSETS_LIST
  do
    UC_OFFSETS=$(echo $OFFSETS_ | tr "[:lower:]" "[:upper:]")
    KOKKOSKERNELS_OFFSETS_CMD="-DKokkosKernels_INST_OFFSET_${UC_OFFSETS}=ON ${KOKKOSKERNELS_OFFSETS_CMD}"
  done
}

get_kernels_layouts_list() {
  echo "parsing KOKKOSKERNELS_LAYOUTS=$KOKKOSKERNELS_LAYOUTS"
  KOKKOSKERNELS_LAYOUTS_CMD=
  PARSE_LAYOUTS_LIST=$(echo $KOKKOSKERNELS_LAYOUTS | tr "," "\n")
  for LAYOUTS_ in $PARSE_LAYOUTS_LIST
  do
    UC_LAYOUTS=$(echo $LAYOUTS_ | tr "[:lower:]" "[:upper:]")
    KOKKOSKERNELS_LAYOUTS_CMD="-DKokkosKernels_INST_${UC_LAYOUTS}=ON ${KOKKOSKERNELS_LAYOUTS_CMD}"
  done
}

get_kernels_spaces_list() {
  echo "parsing KOKKOSKERNELS_SPACES=$KOKKOSKERNELS_SPACES"
  KOKKOSKERNELS_SPACES_CMD=
  PARSE_SPACES_LIST=$(echo $KOKKOSKERNELS_SPACES | tr "," "\n")
  for SPACES_ in $PARSE_SPACES_LIST
  do
    UC_SPACES=$(echo $SPACES_ | tr "[:lower:]" "[:upper:]")
    KOKKOSKERNELS_SPACES_CMD="-DKokkosKernels_INST_MEMSPACE_${UC_SPACES}=ON ${KOKKOSKERNELS_SPACES_CMD}"
  done
}

get_kernels_tpls_list() {
  echo "parsing KOKKOSKERNELS_TPLS=$KOKKOSKERNELS_TPLS"
  KOKKOSKERNELS_TPLS_LIST_CMD=
  KOKKOSKERNELS_USER_TPL_PATH_CMD=
  KOKKOSKERNELS_USER_TPL_LIBNAME_CMD=
  CUBLAS_DEFAULT=OFF
  CUSPARSE_DEFAULT=OFF
  PARSE_TPLS_LIST=$(echo $KOKKOSKERNELS_TPLS | tr "," "\n")
  for TPLS_ in $PARSE_TPLS_LIST
  do
    UC_TPLS=$(echo $TPLS_ | tr "[:lower:]" "[:upper:]")
    KOKKOSKERNELS_TPLS_CMD="-DKokkosKernels_ENABLE_TPL_${UC_TPLS}=ON ${KOKKOSKERNELS_TPLS_CMD}"
    if [ "$UC_TPLS" == "CUBLAS" ]; then
      CUBLAS_DEFAULT=ON
    fi
    if [ "$UC_TPLS" == "CUSPARSE" ]; then
      CUSPARSE_DEFAULT=ON
    fi
    if [ "$UC_TPLS" == "BLAS" ]; then
      if [ "$BLAS_PATH" != "" ]; then
        echo User BLAS_PATH=$BLAS_PATH
        KOKKOSKERNELS_USER_TPL_PATH_CMD="-DBLAS_LIBRARY_DIRS=${BLAS_PATH} ${KOKKOSKERNELS_USER_TPL_PATH_CMD}"
      fi
      if [ "$BLAS_LIBNAME" != "" ]; then
        echo User BLAS_LIBNAME=$BLAS_LIBNAME
        KOKKOSKERNELS_USER_TPL_LIBNAME_CMD="-DBLAS_LIBRARIES=${BLAS_LIBNAME} ${KOKKOSKERNELS_USER_TPL_LIBNAME_CMD}"
      fi
    fi
    if [ "$UC_TPLS" == "LAPACK" ] || [ "$UC_TPLS" == "BLAS" ]; then
      if [ "$LAPACK_PATH" != "" ]; then
        echo User LAPACK_PATH=$LAPACK_PATH
        KOKKOSKERNELS_USER_TPL_PATH_CMD="-DLAPACK_LIBRARY_DIRS=${LAPACK_PATH} ${KOKKOSKERNELS_USER_TPL_PATH_CMD}"
      fi
      if [ "$LAPACK_LIBNAME" != "" ]; then
        echo User LAPACK_LIBNAME=$LAPACK_LIBNAME
        KOKKOSKERNELS_USER_TPL_LIBNAME_CMD="-DLAPACK_LIBRARIES=${LAPACK_LIBNAME} ${KOKKOSKERNELS_USER_TPL_LIBNAME_CMD}"
      fi
    fi
  done
  if [ "$CUBLAS_DEFAULT" == "OFF" ]; then
    KOKKOSKERNELS_TPLS_CMD="-DKokkosKernels_ENABLE_TPL_CUBLAS=OFF ${KOKKOSKERNELS_TPLS_CMD}"
  fi
  if [ "$CUSPARSE_DEFAULT" == "OFF" ]; then
    KOKKOSKERNELS_TPLS_CMD="-DKokkosKernels_ENABLE_TPL_CUSPARSE=OFF ${KOKKOSKERNELS_TPLS_CMD}"
  fi
}

get_kernels_extra_linker_flags() {
  echo "parsing KOKKOSKERNELS_EXTRA_LINKER_FLAGS=$KOKKOSKERNELS_EXTRA_LINKER_FLAGS"
  KOKKOSKERNELS_EXTRA_LINKER_FLAGS_CMD=
  PARSE_EXTRA_LINKER_FLAGS_LIST=$(echo $KOKKOSKERNELS_EXTRA_LINKER_FLAGS | tr "," " ")
  KOKKOSKERNELS_EXTRA_LINKER_FLAGS_CMD="-DCMAKE_EXE_LINKER_FLAGS=${PARSE_EXTRA_LINKER_FLAGS_LIST} ${KOKKOSKERNELS_EXTRA_LINKER_FLAGS_CMD}"
}

display_help_text() {

      echo "KokkosKernels and Kokkos configure options:"
      echo ""
      echo "--kokkos-path=/Path/To/Kokkos:                Path to the Kokkos root directory."
      echo "--kokkos-prefix=/Install/PathToKokkos:        Path to the Kokkos install directory."
      echo "--kokkoskernels-path=/Path/To/KokkosKernels:  Path to the KokkosKernels root directory."
      echo "--prefix=/Install/Path:                       Path to install the KokkosKernels library."
      echo ""
      echo "--with-cuda[=/Path/To/Cuda]:                  Enable Cuda and set path to Cuda Toolkit."
      echo "--with-openmp:                                Enable OpenMP backend."
      echo "--with-pthread:                               Enable Pthreads backend."
      echo "--with-serial:                                Enable Serial backend."
      echo "--with-devices:                               Explicitly add a set of backends."
      echo ""
      echo "--arch=[OPT]:  Set target architectures. Options are:"
      echo "               [AMD]"
      echo "                 AMDAVX          = AMD CPU"
      echo "                 EPYC            = AMD EPYC Zen-Core CPU"
      echo "               [ARM]"
      echo "                 ARMv80          = ARMv8.0 Compatible CPU"
      echo "                 ARMv81          = ARMv8.1 Compatible CPU"
      echo "                 ARMv8-ThunderX  = ARMv8 Cavium ThunderX CPU"
      echo "                 ARMv8-TX2       = ARMv8 Cavium ThunderX2 CPU"
      echo "               [IBM]"
      echo "                 BGQ             = IBM Blue Gene Q"
      echo "                 Power7          = IBM POWER7 and POWER7+ CPUs"
      echo "                 Power8          = IBM POWER8 CPUs"
      echo "                 Power9          = IBM POWER9 CPUs"
      echo "               [Intel]"
      echo "                 WSM             = Intel Westmere CPUs"
      echo "                 SNB             = Intel Sandy/Ivy Bridge CPUs"
      echo "                 HSW             = Intel Haswell CPUs"
      echo "                 BDW             = Intel Broadwell Xeon E-class CPUs"
      echo "                 SKX             = Intel Sky Lake Xeon E-class HPC CPUs (AVX512)"
      echo "               [Intel Xeon Phi]"
      echo "                 KNC             = Intel Knights Corner Xeon Phi"
      echo "                 KNL             = Intel Knights Landing Xeon Phi"
      echo "               [NVIDIA]"
      echo "                 Kepler30        = NVIDIA Kepler generation CC 3.0"
      echo "                 Kepler32        = NVIDIA Kepler generation CC 3.2"
      echo "                 Kepler35        = NVIDIA Kepler generation CC 3.5"
      echo "                 Kepler37        = NVIDIA Kepler generation CC 3.7"
      echo "                 Maxwell50       = NVIDIA Maxwell generation CC 5.0"
      echo "                 Maxwell52       = NVIDIA Maxwell generation CC 5.2"
      echo "                 Maxwell53       = NVIDIA Maxwell generation CC 5.3"
      echo "                 Pascal60        = NVIDIA Pascal generation CC 6.0"
      echo "                 Pascal61        = NVIDIA Pascal generation CC 6.1"
      echo "                 Volta70         = NVIDIA Volta generation CC 7.0"
      echo "                 Volta72         = NVIDIA Volta generation CC 7.2"
      echo ""
      echo "--compiler=/Path/To/Compiler  Set the compiler."
      echo ""
      echo "--debug,-dbg:                 Enable KokkosKernels Debugging."
      echo "--kokkos-debug,-kdbg:          Enable Kokkos Debugging."
      echo "--release:                    Enable KokkosKernels Release Mode."
      echo "--kokkos-release:             Enable Kokkos Release Mode."
      echo ""
      echo "--cxxflags=[FLAGS]            Overwrite CXXFLAGS for library build and test"
      echo "                                build.  This will still set certain required"
      echo "                                flags via KOKKOS_CXXFLAGS (such as -fopenmp,"
      echo "                                --std=c++11, etc.)."
      echo "--cxxstandard=[FLAGS]         Overwrite KOKKOS_CXX_STANDARD for library build and test"
      echo "                                c++11 (default), c++14, c++17, c++1y, c++1z, c++2a"
      echo "--ldflags=[FLAGS]             Overwrite LDFLAGS for library build and test"
      echo "                                build. This will still set certain required"
      echo "                                flags via KOKKOS_LDFLAGS (such as -fopenmp,"
      echo "                                -lpthread, etc.)."
      echo "--with-gtest=/Path/To/Gtest:  Set path to gtest.  (Used in unit and performance"
      echo "                                tests.)"
      echo "--with-hwloc=/Path/To/Hwloc:  Set path to hwloc library."
      echo "--with-memkind=/Path/To/MemKind:  Set path to memkind library."
      echo "--with-options=[OPT]:         Additional options to Kokkos:"
      echo "                                compiler_warnings"
      echo "                                aggressive_vectorization = add ivdep on loops"
      echo "                                disable_profiling = do not compile with profiling hooks"
      echo "                                "
      echo "--with-cuda-options=[OPT]:    Additional options to CUDA:"
      echo "                                force_uvm, use_ldg, enable_lambda, rdc"
      echo "--with-scalars=[SCALARS]:     Set scalars to be instantiated."
      echo "                                Options: float, double, complex_float, complex_double"
      echo "--with-ordinals=[ORDINALS]:   Set ordinals to be instantiated."
      echo "                                Options: int, int64_t"
      echo "--with-offsets=[OFFSETS]:     Set offsets to be instantiated."
      echo "                                Options: int, size_t"
      echo "--with-layouts=[LAYOUTS]:     Set layouts to be instantiated."
      echo "                                Options: layoutleft, layoutright"
      echo "--with-spaces=[SPACES]:       Set spaces to be instantiated."
      echo "                                Options: hostspace, cudaspace, cudauvmspace"
      echo "--with-tpls=[TPLS]:           Set tpls to be instantiated (Proper support requies that appropriate compiler and device must be enabled)."
      echo "                              This may require providing paths and the library name if using custom installs not on a default path"
      echo "                              that CMake searches"
      echo "                                Options: blas, mkl, cublas, cusparse, magma"
      echo "--user-blas-path=[PATH]:      Set path to location of user-specified BLAS library."
      echo "--user-blas-lib=[LIB]:        Library name of desired BLAS install."
      echo "                                Example: For the typical \"libblas.a\" provide \"blas\""
      echo "--user-lapack-path=[PATH]:    Set path to location of user-specified BLAS library."
      echo "--user-lapack-lib=[LIB]:      Library name of the desired LAPACK install."
      echo "                                Example: For the typical \"liblapack.a\" provide \"lapack\""
      echo "--extra-linker-flags=[FLAG]:  Extra linker flags to pass to CMAKE_EXE_LINKER_FLAGS."
      echo "                                Pass flags, separate by comma"
      echo "                                Example: \"-lgfortran,-lma,-Wall\""
#      echo "--with-hpx-options=[OPT]:     Additional options to HPX:"
#      echo "                                enable_async_dispatch"
      echo "--gcc-toolchain=/Path/To/GccRoot:  Set the gcc toolchain to use with clang (e.g. /usr)" 
      echo "--kokkos-make-j=[NUM]:        Set -j parallel level for kokkos install"
      echo "                                Default: j == 4"

}

KOKKOS_INSTALL_PATH=""

KOKKOS_DO_TESTS=OFF
KOKKOS_DO_EXAMPLES=OFF
KOKKOSKERNELS_DO_TESTS=ON

KOKKOS_MAKEINSTALL_J=4

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
    --prefix*)
      PREFIX="${key#*=}"
      ;;
    --kokkos-prefix*)
      KOKKOS_INSTALL_PATH="${key#*=}"
      ;;
    --hpx-path*)
      HPX_PATH="${key#*=}"
      ;;
    --with-cuda)
      update_kokkos_devices Cuda
      CUDA_PATH_NVCC=$(command -v nvcc)
      CUDA_PATH=${CUDA_PATH_NVCC%/bin/nvcc}
      ;;
    # Catch this before '--with-cuda*'
    --with-cuda-options*)
      KOKKOS_CUDA_OPTIONS="${key#*=}"
      ;;
    --with-cuda*)
      update_kokkos_devices Cuda
      CUDA_PATH="${key#*=}"
      ;;
    --with-openmp)
      update_kokkos_devices OpenMP
      ;;
    --with-pthread)
      update_kokkos_devices Pthread
      ;;
    --with-serial)
      update_kokkos_devices Serial
      ;;
    --with-hpx-options*)
      KOKKOS_HPX_OPT="${key#*=}"
      ;;
    --with-hpx*)
      update_kokkos_devices HPX
      if [ -z "$HPX_PATH" ]; then
        HPX_PATH="${key#*=}"
      fi
      ;;
    --with-devices*)
      DEVICES="${key#*=}"
      PARSE_DEVICES=$(echo $DEVICES | tr "," "\n")
      for DEVICE_ in $PARSE_DEVICES
      do 
         update_kokkos_devices $DEVICE_
      done
      ;;
    --with-gtest*)
      GTEST_PATH="${key#*=}"
      ;;
    --with-hwloc*)
      HWLOC_PATH="${key#*=}"
      ;;
    --with-memkind*)
      MEMKIND_PATH="${key#*=}"
      ;;
    --arch*)
      KOKKOS_ARCH="${key#*=}"
      ;;
    --cxxflags*)
      KOKKOS_CXXFLAGS="${key#*=}"
      KOKKOS_CXXFLAGS=${KOKKOS_CXXFLAGS//,/ }
      ;;
    --cxxstandard*)
      KOKKOS_CXX_STANDARD="${key#*=}"
      ;;
    --ldflags*)
      KOKKOS_LDFLAGS="${key#*=}"
      ;;
    --kokkos-debug|-kdbg)
      KOKKOS_DEBUG=yes
      ;;
    --debug|-dbg)
      KOKKOSKERNELS_DEBUG=yes
      ;;
    --kokkos-release)
      KOKKOS_RELEASE=yes
      ;;
    --release)
      KOKKOSKERNELS_RELEASE=yes
      ;;
    --kokkos-make-j*)
      echo "${key} parallel level for kokkos install"
      KOKKOS_MAKEINSTALL_J="${key#*=}"
      ;;
    --enable-kokkos-tests)
      KOKKOS_DO_TESTS=ON
      ;;
    --disable-kokkos-tests)
      # This is the default
      KOKKOS_DO_TESTS=OFF
      ;;
    --enable-kokkos-examples)
      KOKKOS_DO_EXAMPLES=ON
      ;;
    --disable-kokkos-examples)
      # This is the default
      KOKKOS_DO_EXAMPLES=OFF
      ;;
    --enable-tests)
      # This is the default
      KOKKOSKERNELS_DO_TESTS=ON
      ;;
    --disable-tests)
      KOKKOSKERNELS_DO_TESTS=OFF
      ;;
    --compiler*)
      COMPILER="${key#*=}"
      CNUM=$(command -v ${COMPILER} 2>&1 >/dev/null | grep "no ${COMPILER}" | wc -l)
      if [ ${CNUM} -gt 0 ]; then
        echo "Invalid compiler by --compiler command: '${COMPILER}'"
        exit
      fi
      if [[ ! -n  ${COMPILER} ]]; then
        echo "Empty compiler specified by --compiler command."
        exit
      fi
      CNUM=$(command -v ${COMPILER} | grep ${COMPILER} | wc -l)
      if [ ${CNUM} -eq 0 ]; then
        echo "Invalid compiler by --compiler command: '${COMPILER}'"
        exit
      fi
      # ... valid compiler, ensure absolute path set 
      WCOMPATH=$(command -v $COMPILER)
      COMPDIR=$(dirname $WCOMPATH)
      COMPNAME=$(basename $WCOMPATH)
      COMPILER=${COMPDIR}/${COMPNAME}
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
    --with-layouts*)
      KOKKOSKERNELS_LAYOUTS="${key#*=}"
      ;;
    --with-spaces*)
      KOKKOSKERNELS_SPACES="${key#*=}"
      ;;
    --with-tpls*)
      KOKKOSKERNELS_TPLS="${key#*=}"
      ;;
    --user-blas-path*)
      BLAS_PATH="${key#*=}"
      ;;
    --user-blas-lib*)
      BLAS_LIBNAME="${key#*=}"
      ;;
    --user-lapack-path*)
      LAPACK_PATH="${key#*=}"
      ;;
    --user-lapack-lib*)
      LAPACK_LIBNAME="${key#*=}"
      ;;
    --extra-linker-flags*)
      KOKKOSKERNELS_EXTRA_LINKER_FLAGS="${key#*=}"
      ;;
    --with-options*)
      KOKKOS_OPTIONS="${key#*=}"
      ;;
    --gcc-toolchain*)
      KOKKOS_GCC_TOOLCHAIN="${key#*=}"
      ;;
    --help)
      display_help_text
      exit 0
      ;;
    *)
      echo "warning: ignoring unknown option $key"
      ;;
  esac

  shift
done

if [ "$KOKKOS_CXX_STANDARD" == "" ]; then
    STANDARD_CMD=
else
    STANDARD_CMD=-DKokkos_CXX_STANDARD=${KOKKOS_CXX_STANDARD}
fi

if [ "$COMPILER" == "" ]; then
    COMPILER_CMD=
else
    COMPILER_CMD=-DCMAKE_CXX_COMPILER=$COMPILER
fi

echo "KOKKOS_DEBUG CHECK"
if [ "$KOKKOS_DEBUG" == "yes" ]; then
    KOKKOS_BUILDTYPE_CMD=-DCMAKE_BUILD_TYPE=DEBUG
elif [ "$KOKKOS_RELEASE" == "yes" ]; then
    KOKKOS_BUILDTYPE_CMD=-DCMAKE_BUILD_TYPE=RELEASE
else
    KOKKOS_BUILDTYPE_CMD=
fi

echo "KOKKOSKERNELS_DEBUG CHECK"
if [ "$KOKKOSKERNELS_DEBUG" == "yes" ]; then
    KOKKOSKERNELS_BUILDTYPE_CMD=-DCMAKE_BUILD_TYPE=DEBUG
elif [ "$KOKKOSKERNELS_RELEASE" == "yes" ]; then
    KOKKOSKERNELS_BUILDTYPE_CMD=-DCMAKE_BUILD_TYPE=RELEASE
else
    KOKKOSKERNELS_BUILDTYPE_CMD=
fi


# Currently assumes script is in base kokkos-kernels directory
if [ ! -e ${KOKKOSKERNELS_PATH}/CMakeLists.txt ]; then
   if [ "${KOKKOSKERNELS_PATH}" == "" ]; then
   echo "CHECKING: $KOKKOSKERNELS_PATH"
      CM_SCRIPT=$0
      KOKKOSKERNELS_PATH=`dirname $CM_SCRIPT`
      if [ ! -e ${KOKKOSKERNELS_PATH}/CMakeLists.txt ]; then
         echo "${KOKKOSKERNELS_PATH} repository appears to not be complete.  please verify and try again"
         exit 0
      fi
   else
      echo "KOKKOSKERNELS_PATH does not appear to be set properly. please specify in location of CMakeLists.txt"   
      display_help_text
      exit 0
   fi
fi

if [ "${KOKKOS_PATH}" == "" ]; then
  CM_SCRIPT=$0
  KOKKOS_PATH=`dirname $CM_SCRIPT`
  KOKKOS_PATH="${KOKKOS_PATH}/../kokkos"
  echo "CHECKING: $KOKKOS_PATH"
  $(ls $KOKKOS_PATH)
  if [ ! -e ${KOKKOS_PATH}/CMakeLists.txt ]; then
     echo "Either kokkos repository is not in the same base directory as kokkos-kernels or ${KOKKOS_PATH} repository appears to not be complete.  Please verify or provide the path to kokkos and try again"
     exit 0
  else
     echo "${KOKKOS_PATH} repository appears complete."
  fi
else
  # IMPROVE THIS CHECK - not sufficient for old vs new kokkos
  if [ ! -e ${KOKKOS_PATH}/CMakeLists.txt ]; then
     echo "KOKKOS_PATH does not appear to be set properly. Please check and try again."   
     display_help_text
     exit 0
  fi
fi

get_kokkos_device_list
get_kokkos_option_list
get_kokkos_arch_list
get_kokkos_cuda_option_list

get_kernels_scalar_list
get_kernels_ordinals_list
get_kernels_offsets_list
get_kernels_layouts_list
get_kernels_spaces_list
get_kernels_tpls_list
get_kernels_extra_linker_flags

## if HPX is enabled, we need to enforce cxx standard = 14
if [[ ${KOKKOS_DEVICE_CMD} == *Kokkos_ENABLE_HPX* ]]; then
   if [ ${#KOKKOS_CXX_STANDARD} -lt 14 ]; then
      echo CXX Standard must be 14 or higher for HPX to work.
      KOKKOS_CXX_STANDARD=14
   fi
fi

if [[ ${COMPILER} == *clang* ]]; then
   gcc_path=$(which g++ | awk --field-separator='/bin/g++' '{printf $1}' )
   KOKKOS_CXXFLAGS="${KOKKOS_CXXFLAGS} --gcc-toolchain=${gcc_path}"

   if [ ! "${CUDA_PATH}" == "" ]; then
      KOKKOS_CXXFLAGS="${KOKKOS_CXXFLAGS} --cuda-path=${CUDA_PATH}"
   fi 
fi


KOKKOS_INSTALL_DIRNAME="kokkos-install"
#if [ "${PREFIX}" == "" ]; then
#fi
if [ "${KOKKOS_INSTALL_PATH}" == "" ]; then
  KOKKOS_INSTALL_PATH="${PWD}/$KOKKOS_INSTALL_DIRNAME"
  echo "KOKKOS_INSTALL_PATH=$KOKKOS_INSTALL_PATH"
else
#  User-provided install path
#  KOKKOS_INSTALL_PATH="${PREFIX}/$KOKKOS_INSTALL_DIRNAME"
  echo "KOKKOS_INSTALL_PATH=$KOKKOS_INSTALL_PATH"
fi

STORE_KOKKOSKERNELS_BUILD_PATH=${PWD}

# Build Kokkos
mkdir -p ${KOKKOS_INSTALL_PATH}
cd ${KOKKOS_INSTALL_PATH}

#KOKKOS_INSTALL_PATH="${PWD}"

echo ""
echo cmake $COMPILER_CMD  -DCMAKE_CXX_FLAGS="${KOKKOS_CXXFLAGS}" -DCMAKE_EXE_LINKER_FLAGS="${KOKKOS_LDFLAGS}" -DCMAKE_INSTALL_PREFIX=${KOKKOS_INSTALL_PATH} ${KOKKOS_DEVICE_CMD} ${KOKKOS_ARCH_CMD} -DKokkos_ENABLE_TESTS=${KOKKOS_DO_TESTS} -DKokkos_ENABLE_EXAMPLES=${KOKKOS_DO_EXAMPLES} ${KOKKOS_OPTION_CMD} ${KOKKOS_CUDA_OPTION_CMD} -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_CXX_EXTENSIONS=OFF ${STANDARD_CMD} ${KOKKOS_BUILDTYPE_CMD} ${KOKKOS_PATH}
echo ""
cmake $COMPILER_CMD  -DCMAKE_CXX_FLAGS="${KOKKOS_CXXFLAGS//\"}" -DCMAKE_EXE_LINKER_FLAGS="${KOKKOS_LDFLAGS//\"}" -DCMAKE_INSTALL_PREFIX=${KOKKOS_INSTALL_PATH} ${KOKKOS_DEVICE_CMD} ${KOKKOS_ARCH_CMD} -DKokkos_ENABLE_TESTS=${KOKKOS_DO_TESTS} -DKokkos_ENABLE_EXAMPLES=${KOKKOS_DO_EXAMPLES} ${KOKKOS_OPTION_CMD} ${KOKKOS_CUDA_OPTION_CMD} -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_CXX_EXTENSIONS=OFF ${STANDARD_CMD} ${KOKKOS_BUILDTYPE_CMD} ${KOKKOS_PATH}

# Install kokkos library
make install -j $KOKKOS_MAKEINSTALL_J

# Replacing Kokkos_DIR with Kokkos_ROOT may resolve need for this check...
KOKKOS_FIND_PATH=
if [ -d ${KOKKOS_INSTALL_PATH}/lib/cmake/Kokkos ]; then
  KOKKOS_FIND_PATH=${KOKKOS_INSTALL_PATH}/lib/cmake/Kokkos
elif [ -d ${KOKKOS_INSTALL_PATH}/lib64/cmake/Kokkos ]; then
  KOKKOS_FIND_PATH=${KOKKOS_INSTALL_PATH}/lib64/cmake/Kokkos
else
  echo "Error: Kokkos install did not yield <kokkos-install-path>/lib/cmake/Kokkos or <kokkos-install-path>/lib64/cmake/Kokkos"
  exit
fi

#echo "KOKKOS_INSTALL_PATH = ${KOKKOS_INSTALL_PATH}"
#echo "KOKKOKS_FIND_PATH = ${KOKKOS_FIND_PATH}"
#echo "STORE_KOKKOSKERNELS_BUILD_PATH  = ${STORE_KOKKOSKERNELS_BUILD_PATH}"
#echo "KOKKOS_INSTALL_PATH: ${KOKKOS_INSTALL_PATH}"

cd $STORE_KOKKOSKERNELS_BUILD_PATH

echo ""
echo cmake $COMPILER_CMD -DKokkos_DIR="${KOKKOS_FIND_PATH}" -DCMAKE_CXX_FLAGS="${KOKKOS_CXXFLAGS}" -DCMAKE_INSTALL_PREFIX="${PREFIX}" -DKokkosKernels_ENABLE_TESTS=${KOKKOSKERNELS_DO_TESTS} ${KOKKOSKERNELS_SCALARS_CMD} ${KOKKOSKERNELS_ORDINALS_CMD} ${KOKKOSKERNELS_OFFSETS_CMD} ${KOKKOSKERNELS_LAYOUTS_CMD} ${KOKKOSKERNELS_TPLS_CMD} ${KOKKOSKERNELS_USER_TPL_PATH_CMD} ${KOKKOSKERNELS_USER_TPL_LIBNAME_CMD} ${KOKKOSKERNELS_EXTRA_LINKER_FLAGS_CMD} ${KOKKOSKERNELS_BUILDTYPE_CMD} ${KOKKOSKERNELS_SPACES_CMD} ${KOKKOSKERNELS_PATH}
echo ""
cmake $COMPILER_CMD -DKokkos_DIR="${KOKKOS_FIND_PATH}" -DCMAKE_CXX_FLAGS="${KOKKOS_CXXFLAGS//\"}" -DCMAKE_INSTALL_PREFIX="${PREFIX}" -DKokkosKernels_ENABLE_TESTS=${KOKKOSKERNELS_DO_TESTS} ${KOKKOSKERNELS_SCALARS_CMD} ${KOKKOSKERNELS_ORDINALS_CMD} ${KOKKOSKERNELS_OFFSETS_CMD} ${KOKKOSKERNELS_LAYOUTS_CMD} ${KOKKOSKERNELS_TPLS_CMD} ${KOKKOSKERNELS_USER_TPL_PATH_CMD} ${KOKKOSKERNELS_USER_TPL_LIBNAME_CMD} ${KOKKOSKERNELS_EXTRA_LINKER_FLAGS_CMD} ${KOKKOSKERNELS_BUILDTYPE_CMD} ${KOKKOSKERNELS_SPACES_CMD} ${KOKKOSKERNELS_PATH}

