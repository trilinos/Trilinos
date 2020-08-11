#!/bin/bash

CURRENT_SCRIPTS_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${CURRENT_SCRIPTS_DIR}/../..`

#
# Test compiler parsing
#

# Make work on all systems reguardless of ATDM Trilinos env loaded
unset ATDM_CONFIG_SYSTEM_NAME
unset ATDM_CONFIG_SYSTEM_DIR


testAllDefaults() {
  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} DEFAULT
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_KOKKOS_ARCH} DEFAULT
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} DEBUG
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_NODE_TYPE} SERIAL
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_OPENMP} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_CUDA} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_PTHREADS} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_FPIC} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPLEX} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_PT_PACKAGES} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_FINISHED_SET_BUILD_OPTIONS} 1
}


testCompilerClangAndDefaults() {
  ATDM_CONFIG_BUILD_NAME=clang
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_KOKKOS_ARCH} DEFAULT
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} DEBUG
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_NODE_TYPE} SERIAL
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_OPENMP} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_CUDA} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_PTHREADS} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_FPIC} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPLEX} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_PT_PACKAGES} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_FINISHED_SET_BUILD_OPTIONS} 1
}


testCompilerClang() {

  ATDM_CONFIG_BUILD_NAME=clang-5.0.1
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-5.0.1

  ATDM_CONFIG_BUILD_NAME=before_clang-5.0.1_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-5.0.1

  ATDM_CONFIG_BUILD_NAME=before-clang-5.0.1-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-5.0.1

  ATDM_CONFIG_BUILD_NAME=clang-7.0.1
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CLANG-7.0.1

}


testCompilerCuda() {

  ATDM_CONFIG_BUILD_NAME=cuda
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CUDA

  ATDM_CONFIG_BUILD_NAME=cuda-9.2
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CUDA-9.2

  ATDM_CONFIG_BUILD_NAME=cuda-9.2-gnu-7.2.0
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CUDA-9.2_GNU-7.2.0

  ATDM_CONFIG_BUILD_NAME=cuda-9.2_gnu-7.2.0
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CUDA-9.2_GNU-7.2.0

  ATDM_CONFIG_BUILD_NAME=cuda-10.0-gnu-7.4.0
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} CUDA-10.0_GNU-7.4.0

}


testCompilerIntel() {

  ATDM_CONFIG_BUILD_NAME=intel
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} INTEL

  ATDM_CONFIG_BUILD_NAME=intel-17
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} INTEL-17.0.1

  ATDM_CONFIG_BUILD_NAME=intel-17.0.1
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} INTEL-17.0.1

  ATDM_CONFIG_BUILD_NAME=intel-18
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} INTEL-18.0.5

}


testNompi() {

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_MPI} ON

  ATDM_CONFIG_BUILD_NAME=default-mpi
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_MPI} ON

  ATDM_CONFIG_BUILD_NAME=default-MPI
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_MPI} ON

  ATDM_CONFIG_BUILD_NAME=default-mpi-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_MPI} ON

  ATDM_CONFIG_BUILD_NAME=default-no-mpi
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_MPI} OFF

  ATDM_CONFIG_BUILD_NAME=default-NO-MPI
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_MPI} OFF

  ATDM_CONFIG_BUILD_NAME=default-no-mpi-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_MPI} OFF

}


testKokkosArch() {

  # Test first arch and defaults for everything else
  ATDM_CONFIG_BUILD_NAME=default-AMDAVX
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPILER} DEFAULT
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_KOKKOS_ARCH} AMDAVX
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} DEBUG
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_NODE_TYPE} SERIAL
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_OPENMP} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_CUDA} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_PTHREADS} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_FPIC} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPLEX} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_PT_PACKAGES} OFF

  # Test arch in the middle
  ATDM_CONFIG_BUILD_NAME=default-HSW
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_KOKKOS_ARCH} HSW

  # Test arch at the end
  ATDM_CONFIG_BUILD_NAME=default-TX2
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_KOKKOS_ARCH} TX2

  # Match the first in the list of arch
  ATDM_CONFIG_BUILD_NAME=default-ARMv8-ThunderX-TX2
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_KOKKOS_ARCH} ARMv8-ThunderX

  # Match the first in the list of arch, even if listed later
  ATDM_CONFIG_BUILD_NAME=default-TX2-ARMv8-ThunderX
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_KOKKOS_ARCH} ARMv8-ThunderX

  # At beginning lower -
  ATDM_CONFIG_BUILD_NAME=hsw-default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_KOKKOS_ARCH} HSW

  # At beginning upper _
  ATDM_CONFIG_BUILD_NAME=KNL_default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_KOKKOS_ARCH} KNL

  # In middle lower -
  ATDM_CONFIG_BUILD_NAME=default-hsw-dummy
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_KOKKOS_ARCH} HSW

  # In middle upper _
  ATDM_CONFIG_BUILD_NAME=default_KNL_dummy
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_KOKKOS_ARCH} KNL

  # At end lower -
  ATDM_CONFIG_BUILD_NAME=default-hsw
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_KOKKOS_ARCH} HSW

  # AT end upper _
  ATDM_CONFIG_BUILD_NAME=default_KNL
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_KOKKOS_ARCH} KNL

}


testBuildType() {

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} DEBUG

  ATDM_CONFIG_BUILD_NAME=default-somethingopt
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} DEBUG

  ATDM_CONFIG_BUILD_NAME=default-releases
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} DEBUG

  ATDM_CONFIG_BUILD_NAME=default-optend
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} DEBUG

  ATDM_CONFIG_BUILD_NAME=default-yesrelease
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} DEBUG

  ATDM_CONFIG_BUILD_NAME=default-release-debug
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} RELEASE-DEBUG

  ATDM_CONFIG_BUILD_NAME=default_release_debug
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} RELEASE-DEBUG

  ATDM_CONFIG_BUILD_NAME=default_opt-dbg
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} RELEASE-DEBUG

  ATDM_CONFIG_BUILD_NAME=default-opt_dbg
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} RELEASE-DEBUG

  ATDM_CONFIG_BUILD_NAME=default-release
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} RELEASE 

  ATDM_CONFIG_BUILD_NAME=default-debug
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} DEBUG

  ATDM_CONFIG_BUILD_NAME=default-dbg
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_BUILD_TYPE} DEBUG

}


testNodeType() {

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_NODE_TYPE} SERIAL
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_CUDA} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_OPENMP} OFF

  ATDM_CONFIG_BUILD_NAME=something-cuda-10.2-AFTER
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_NODE_TYPE} CUDA
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_CUDA} ON
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_OPENMP} OFF

  ATDM_CONFIG_BUILD_NAME=default-serial-AFTER
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_NODE_TYPE} SERIAL
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_CUDA} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_OPENMP} OFF

  ATDM_CONFIG_BUILD_NAME=default-openmp-AFTER
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_NODE_TYPE} OPENMP
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_CUDA} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_OPENMP} ON

  ATDM_CONFIG_BUILD_NAME=default-openmpi-AFTER
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_NODE_TYPE} SERIAL
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_CUDA} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_OPENMP} OFF

  ATDM_CONFIG_BUILD_NAME=default-openmpi-4.0.2-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_NODE_TYPE} SERIAL
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_CUDA} OFF
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_USE_OPENMP} OFF

}


testRDC() {

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} OFF

  ATDM_CONFIG_BUILD_NAME=default-rdc
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} ON

  ATDM_CONFIG_BUILD_NAME=default_rdc
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} ON

  ATDM_CONFIG_BUILD_NAME=default-rdcs
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} OFF

  ATDM_CONFIG_BUILD_NAME=default_rdc-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} ON

  ATDM_CONFIG_BUILD_NAME=default_RDC-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} ON

  ATDM_CONFIG_BUILD_NAME=default_rdc_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} ON

  ATDM_CONFIG_BUILD_NAME=default-no-rdc
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} OFF

  ATDM_CONFIG_BUILD_NAME=default_no-rdc
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} OFF

  ATDM_CONFIG_BUILD_NAME=default_rdc_no-rdc
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} OFF

  ATDM_CONFIG_BUILD_NAME=default-rdc-no-rdc
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} OFF

  ATDM_CONFIG_BUILD_NAME=default_rdc_no-rdc-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_CUDA_RDC} OFF

}


testFPIC() {

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_FPIC} OFF

  ATDM_CONFIG_BUILD_NAME=default-fpic-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_FPIC} ON

  ATDM_CONFIG_BUILD_NAME=default_fpic-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_FPIC} ON

}


testComplex() {

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPLEX} OFF

  ATDM_CONFIG_BUILD_NAME=default-complex
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPLEX} ON

  ATDM_CONFIG_BUILD_NAME=default-COMPLEX
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPLEX} ON

  ATDM_CONFIG_BUILD_NAME=default-complexs
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPLEX} OFF

  ATDM_CONFIG_BUILD_NAME=default_complex_after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPLEX} ON

  ATDM_CONFIG_BUILD_NAME=default-no-complex
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPLEX} OFF

  ATDM_CONFIG_BUILD_NAME=default-nos-complex
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPLEX} ON

  ATDM_CONFIG_BUILD_NAME=default-no-complex-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPLEX} OFF

  ATDM_CONFIG_BUILD_NAME=default-no-COMPLEX-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_COMPLEX} OFF

}


testSharedLibs() {

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} OFF

  ATDM_CONFIG_BUILD_NAME=default-static
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} OFF

  ATDM_CONFIG_BUILD_NAME=default-static-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} OFF

  ATDM_CONFIG_BUILD_NAME=default-shared
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} ON

  ATDM_CONFIG_BUILD_NAME=default-shared-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} ON

  ATDM_CONFIG_BUILD_NAME=default-static-shared-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} ON

  ATDM_CONFIG_BUILD_NAME=default-shared-static-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} ON

  ATDM_CONFIG_BUILD_NAME=default-static-SHARED-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} ON

  ATDM_CONFIG_BUILD_NAME=default-static-SHAREDs-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_SHARED_LIBS} OFF

}


testPrimaryTested() {

  ATDM_CONFIG_BUILD_NAME=default
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_PT_PACKAGES} OFF

  ATDM_CONFIG_BUILD_NAME=default-pt
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_PT_PACKAGES} ON

  ATDM_CONFIG_BUILD_NAME=default-PT
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_PT_PACKAGES} ON

  ATDM_CONFIG_BUILD_NAME=default-pt-after
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_PT_PACKAGES} ON

  ATDM_CONFIG_BUILD_NAME=default-pts
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_PT_PACKAGES} OFF

  ATDM_CONFIG_BUILD_NAME=default_pt
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_PT_PACKAGES} ON

  ATDM_CONFIG_BUILD_NAME=default_PT-fater
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_PT_PACKAGES} ON

  ATDM_CONFIG_BUILD_NAME=default_PTS-fater
  . ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh
  ${_ASSERT_EQUALS_} ${ATDM_CONFIG_PT_PACKAGES} OFF

}


#
# Run the unit tests
#

SHUNIT2_DIR=`readlink -f ${ATDM_CONFIG_SCRIPT_DIR}/../../../commonTools/test/shunit2`
. ${SHUNIT2_DIR}/shunit2
