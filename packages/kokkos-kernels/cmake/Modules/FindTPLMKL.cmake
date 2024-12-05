find_package(MKL)
IF(TARGET MKL::MKL)
  # MKL version >= 2021 (see kokkos wiki and intel documentation. MKL CMake module file has been introduced starting MKL >= 2021)
  IF (KOKKOS_ENABLE_SYCL) #get from kokkos-core
    # MKL version >= 2022 (see kokkos wiki)
    IF (NOT TARGET MKL::MKL_DPCPP)
      MESSAGE(FATAL_ERROR "KOKKOS_ENABLE_SYCL activated but the target MKL_DPCPP wasn't found")
    ENDIF()
  ENDIF()
  SET(TPL_MKL_IMPORTED_NAME MKL::MKL)
  SET(TPL_IMPORTED_NAME MKL::MKL)
  ADD_LIBRARY(MKL INTERFACE)
  IF(KOKKOS_ENABLE_SYCL)
    TARGET_LINK_LIBRARIES(MKL INTERFACE MKL::MKL MKL::MKL_DPCPP)
  ELSE()
    TARGET_LINK_LIBRARIES(MKL INTERFACE MKL::MKL )
  ENDIF()
  ADD_LIBRARY(KokkosKernels::MKL ALIAS MKL )
  GET_TARGET_PROPERTY(LIB_TYPE ${TPL_IMPORTED_NAME} TYPE)
  MESSAGE("LIB_TYPE: ${LIB_TYPE}")
  # kokkoskernels_export_imported_tpl install MKL with target name MKL instead of
  # MKL::MKL or KokkosKernels::MKL, so we need to install a specific ALIAS one
  if(TARGET MKL)
    MESSAGE("TARGET MKL CREATED")
  ENDIF()
ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
# Regular way with MKL version < 2021 (Where MKL doesn't provide cmake module file)
  TRY_COMPILE(KOKKOSKERNELS_HAS_MKL_ARG
    ${KOKKOSKERNELS_TOP_BUILD_DIR}/tpl_tests
    ${KOKKOSKERNELS_TOP_SOURCE_DIR}/cmake/compile_tests/mkl.cpp
    LINK_LIBRARIES -mkl
    COMPILE_DEFINITIONS -mkl)
  KOKKOSKERNELS_CREATE_IMPORTED_TPL(MKL INTERFACE COMPILE_OPTIONS -mkl LINK_OPTIONS -mkl)
  INCLUDE(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(TPLMKL DEFAULT_MSG KOKKOSKERNELS_HAS_MKL_ARG)
ELSEIF(WIN32)
  SET(BLA_VENDOR Intel10_64lp)
  FIND_PACKAGE(BLAS REQUIRED)
  IF (NOT DEFINED ENV{MKLROOT})
    SET(NO_MKL_ROOT_GIVEN "MKL-NOTFOUND")
    MESSAGE(WARNING "No MKLROOT environment variable specified - must source mklvars.sh to configure MKL path")
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL NO_MKL_ROOT_GIVEN)
  ELSE()
    KOKKOSKERNELS_CREATE_IMPORTED_TPL(MKL INTERFACE
      LINK_OPTIONS ${BLAS_LINKER_FLAGS}
      LINK_LIBRARIES ${BLAS_LIBRARIES}
    )
  ENDIF()
ELSE()
  IF (NOT DEFINED ENV{MKLROOT})
    SET(NO_MKL_ROOT_GIVEN "MKL-NOTFOUND")
    MESSAGE(WARNING "No MKLROOT environment variable specified - must source mklvars.sh to configure MKL path")
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL NO_MKL_ROOT_GIVEN)
  ELSE()
    SET(MKL_ROOT $ENV{MKLROOT})
    #go ahead and use LD_LIBRARY_PATH to find certain libs
    LIST(APPEND ENV_LIBDIRS ENV LD_LIBRARY_PATH)
    #override what CMake looks for
    #gnu_thread does not work on some platforms
    #just always use intel_thread
    KOKKOSKERNELS_FIND_IMPORTED(MKL INTERFACE
      LIBRARIES
        mkl_intel_lp64
        mkl_intel_thread
        mkl_core
        iomp5
      LIBRARY_PATHS
        ${MKL_ROOT}/lib/intel64
        ${ENV_LIBDIRS}
      HEADER
        mkl.h
      HEADER_PATHS
        ${MKL_ROOT}/include
    )
    ENDIF()
ENDIF()
# This logic to find MKL is only used in non-Trilinos builds.
# In this case, MKL can always be used as the host BLAS/LAPACK implementation
# (whether MKL_INT is 32- or 64-bit).
set (MKL_PROVIDES_BLAS_LAPACK ON INTERNAL)
