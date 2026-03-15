find_package(MKL)
if(TARGET MKL::MKL)
  # MKL version >= 2021 (see kokkos wiki and intel documentation. MKL CMake module file has been introduced starting MKL >= 2021)
  if(KOKKOS_ENABLE_SYCL) #get from kokkos-core
    # MKL version >= 2022 (see kokkos wiki)
    if(NOT TARGET MKL::MKL_DPCPP)
      message(FATAL_ERROR "KOKKOS_ENABLE_SYCL activated but the target MKL_DPCPP wasn't found")
    endif()
  endif()
  set(TPL_MKL_IMPORTED_NAME MKL::MKL)
  set(TPL_IMPORTED_NAME MKL::MKL)
  add_library(MKL INTERFACE)
  if(KOKKOS_ENABLE_SYCL)
    target_link_libraries(MKL INTERFACE MKL::MKL MKL::MKL_DPCPP)
  else()
    target_link_libraries(MKL INTERFACE MKL::MKL)
  endif()
  add_library(KokkosKernels::MKL ALIAS MKL)
  get_target_property(LIB_TYPE ${TPL_IMPORTED_NAME} TYPE)
  message("LIB_TYPE: ${LIB_TYPE}")
  # kokkoskernels_export_imported_tpl install MKL with target name MKL instead of
  # MKL::MKL or KokkosKernels::MKL, so we need to install a specific ALIAS one
  if(TARGET MKL)
    message("TARGET MKL CREATED")
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  # Regular way with MKL version < 2021 (Where MKL doesn't provide cmake module file)
  try_compile(
    KOKKOSKERNELS_HAS_MKL_ARG ${KOKKOSKERNELS_TOP_BUILD_DIR}/tpl_tests
    ${KOKKOSKERNELS_TOP_SOURCE_DIR}/cmake/compile_tests/mkl.cpp
    LINK_LIBRARIES -mkl
    COMPILE_DEFINITIONS -mkl)
  kokkoskernels_create_imported_tpl(MKL INTERFACE COMPILE_OPTIONS -mkl LINK_OPTIONS -mkl)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(TPLMKL DEFAULT_MSG KOKKOSKERNELS_HAS_MKL_ARG)
elseif(WIN32)
  set(BLA_VENDOR Intel10_64lp)
  find_package(BLAS REQUIRED)
  if(NOT DEFINED ENV{MKLROOT})
    set(NO_MKL_ROOT_GIVEN "MKL-NOTFOUND")
    message(WARNING "No MKLROOT environment variable specified - must source mklvars.sh to configure MKL path")
    find_package_handle_standard_args(MKL NO_MKL_ROOT_GIVEN)
  else()
    kokkoskernels_create_imported_tpl(MKL INTERFACE LINK_OPTIONS ${BLAS_LINKER_FLAGS} LINK_LIBRARIES ${BLAS_LIBRARIES})
  endif()
else()
  if(NOT DEFINED ENV{MKLROOT})
    set(NO_MKL_ROOT_GIVEN "MKL-NOTFOUND")
    message(WARNING "No MKLROOT environment variable specified - must source mklvars.sh to configure MKL path")
    find_package_handle_standard_args(MKL NO_MKL_ROOT_GIVEN)
  else()
    set(MKL_ROOT $ENV{MKLROOT})
    #go ahead and use LD_LIBRARY_PATH to find certain libs
    list(APPEND ENV_LIBDIRS ENV LD_LIBRARY_PATH)
    #override what CMake looks for
    #gnu_thread does not work on some platforms
    #just always use intel_thread
    kokkoskernels_find_imported(MKL INTERFACE
      LIBRARIES     mkl_intel_lp64 mkl_intel_thread mkl_core iomp5
      LIBRARY_PATHS ${MKL_ROOT}/lib/intel64 ${ENV_LIBDIRS}
      HEADER        mkl.h
      HEADER_PATHS  ${MKL_ROOT}/include)
  endif()
endif()
# This logic to find MKL is only used in non-Trilinos builds.
# In this case, MKL can always be used as the host BLAS/LAPACK implementation
# (whether MKL_INT is 32- or 64-bit).
set(MKL_PROVIDES_BLAS_LAPACK ON INTERNAL)
