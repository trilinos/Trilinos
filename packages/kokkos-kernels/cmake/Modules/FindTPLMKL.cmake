IF (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  TRY_COMPILE(KOKKOSKERNELS_HAS_MKL_ARG
    ${KOKKOSKERNELS_TOP_BUILD_DIR}/tpl_tests
    ${KOKKOSKERNELS_TOP_SOURCE_DIR}/cmake/compile_tests/mkl.cpp
    LINK_LIBRARIES -mkl
    COMPILE_DEFINITIONS -mkl)
  KOKKOSKERNELS_CREATE_IMPORTED_TPL(MKL INTERFACE COMPILE_OPTIONS -mkl LINK_OPTIONS -mkl)
  INCLUDE(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG KOKKOSKERNELS_HAS_MKL_ARG)
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
