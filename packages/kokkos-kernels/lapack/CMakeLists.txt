list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lapack/src)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lapack/impl)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lapack/eti)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/lapack/eti)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lapack/tpls)

# Adding unit-tests
kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR}/lapack)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/lapack)

#########################
#                       #
# Logic for LAPACK TPLs #
#                       #
#########################

#Include LAPACK, Lapack host wrapper
if(KOKKOSKERNELS_ENABLE_TPL_LAPACK OR KOKKOSKERNELS_ENABLE_TPL_MKL OR KOKKOSKERNELS_ENABLE_TPL_ARMPL OR KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  #Do NOT add this to include path
  append_glob(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tpls/KokkosLapack_Host_tpl.cpp)
endif()

# Include host lapack TPL source file
if(KOKKOSKERNELS_ENABLE_TPL_LAPACK OR KOKKOSKERNELS_ENABLE_TPL_MKL OR KOKKOSKERNELS_ENABLE_TPL_ARMPL OR KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  list(APPEND SOURCES lapack/tpls/KokkosLapack_Host_tpl.cpp)
endif()

# Include cuda lapack TPL source file
if(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)
  list(APPEND SOURCES lapack/tpls/KokkosLapack_Cuda_tpl.cpp)
endif()

if(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
  list(APPEND SOURCES lapack/tpls/KokkosLapack_Magma_tpl.cpp)
endif()

##################
#                #
# ETI generation #
#                #
##################

#Build up a list of DECL, AVAIL, and INST macros
#that should be instantiated based on input options
#Generate @X@ variables in the template X.hpp.in and X.cpp.in
#files containing the list of all needed macros

kokkoskernels_generate_eti(Lapack_gesv gesv
  COMPONENTS  lapack
  HEADER_LIST ETI_HEADERS
  SOURCE_LIST SOURCES
  TYPE_LISTS  FLOATS LAYOUTS DEVICES)

kokkoskernels_generate_eti(Lapack_trtri trtri
  COMPONENTS  lapack
  HEADER_LIST ETI_HEADERS
  SOURCE_LIST SOURCES
  TYPE_LISTS  FLOATS LAYOUTS DEVICES)

kokkoskernels_generate_eti(Lapack_svd svd
  COMPONENTS  lapack
  HEADER_LIST ETI_HEADERS
  SOURCE_LIST SOURCES
  TYPE_LISTS  FLOATS LAYOUTS DEVICES)

KOKKOSKERNELS_GENERATE_ETI(Lapack_geqrf geqrf
  COMPONENTS  lapack
  HEADER_LIST ETI_HEADERS
  SOURCE_LIST SOURCES
  TYPE_LISTS  FLOATS LAYOUTS DEVICES)
