LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lapack/src)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lapack/impl)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lapack/eti)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/lapack/eti)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lapack/tpls)

# Adding unit-tests
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/lapack)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/lapack)

#########################
#                       #
# Logic for LAPACK TPLs #
#                       #
#########################

#Include LAPACK, Lapack host wrapper
IF (KOKKOSKERNELS_ENABLE_TPL_LAPACK OR KOKKOSKERNELS_ENABLE_TPL_MKL OR KOKKOSKERNELS_ENABLE_TPL_ARMPL)
  #Do NOT add this to include path
  APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tpls/KokkosLapack_Host_tpl.cpp)
ENDIF()

# Include host lapack TPL source file
IF (KOKKOSKERNELS_ENABLE_TPL_LAPACK OR KOKKOSKERNELS_ENABLE_TPL_MKL OR KOKKOSKERNELS_ENABLE_TPL_ARMPL)
  LIST(APPEND SOURCES
   lapack/tpls/KokkosLapack_Host_tpl.cpp
  )
ENDIF()

# Include cuda lapack TPL source file
IF (KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)
  LIST(APPEND SOURCES
   lapack/tpls/KokkosLapack_Cuda_tpl.cpp
  )
ENDIF()

IF (KOKKOSKERNELS_ENABLE_TPL_MAGMA)
  LIST(APPEND SOURCES
   lapack/tpls/KokkosLapack_Magma_tpl.cpp
  )
ENDIF()

##################
#                #
# ETI generation #
#                #
##################

#Build up a list of DECL, AVAIL, and INST macros
#that should be instantiated based on input options
#Generate @X@ variables in the template X.hpp.in and X.cpp.in
#files containing the list of all needed macros

KOKKOSKERNELS_GENERATE_ETI(Lapack_gesv gesv
  COMPONENTS  lapack
  HEADER_LIST ETI_HEADERS
  SOURCE_LIST SOURCES
  TYPE_LISTS  FLOATS LAYOUTS DEVICES
)

KOKKOSKERNELS_GENERATE_ETI(Lapack_trtri trtri
  COMPONENTS  lapack
  HEADER_LIST ETI_HEADERS
  SOURCE_LIST SOURCES
  TYPE_LISTS  FLOATS LAYOUTS DEVICES
)

KOKKOSKERNELS_GENERATE_ETI(Lapack_svd svd
  COMPONENTS  lapack
  HEADER_LIST ETI_HEADERS
  SOURCE_LIST SOURCES
  TYPE_LISTS  FLOATS LAYOUTS DEVICES
)
