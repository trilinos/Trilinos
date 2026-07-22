# Adding source directory to the build
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/batched/eti)

list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/src)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/impl)

list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/sparse/src)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/sparse/impl)

if(NOT KokkosKernels_ENABLE_COMPONENT_BLAS)
  list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/blas/impl)
  list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/blas/tpls)
  append_glob(HEADERS ${PACKAGE_SOURCE_DIR}/blas/impl/KokkosBlas_util.hpp)
  list(APPEND SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/batched/KokkosBatched_Util.cpp)
endif()

if(KokkosKernels_ENABLE_TESTS OR KokkosKernels_ENABLE_TESTS_AND_PERFSUITE)
  # Adding unit-tests
  kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR}/batched)
  kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/batched)
endif()
# NOTE: Above, the build directory 'batched' is not created unless unit tests
# are actually enabled (which are actually included from the base-level
# CMakeLists.txt file).  And the KokkosKernelsTargets.cmake file that gets
# generated from this CMake package in the build dir will be broken if these
# are listed in the `INTERFACE_INCLUDE_DIRECTORIES` property when the build
# `batched` is not created (see Trilinos PR #11966).

macro(gen_batched_eti FUNC_NAME)
  kokkoskernels_generate_eti("${FUNC_NAME}_bll" Gemm
    COMPONENTS  batched
    HEADER_LIST ETI_HEADERS
    SOURCE_LIST SOURCES
    TYPE_LISTS  FLOATS RIGHT_LAYOUTS DEVICES)
  kokkoskernels_generate_eti("${FUNC_NAME}_blr" Gemm
    COMPONENTS  batched
    HEADER_LIST ETI_HEADERS
    SOURCE_LIST SOURCES
    TYPE_LISTS  FLOATS LEFT_LAYOUTS DEVICES)
endmacro()

gen_batched_eti(Batched_Gemm_nt_nt)
gen_batched_eti(Batched_Gemm_nt_t)
gen_batched_eti(Batched_Gemm_t_nt)
gen_batched_eti(Batched_Gemm_t_t)
