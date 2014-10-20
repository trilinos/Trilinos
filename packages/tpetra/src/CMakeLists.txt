INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Add Kokkos refactor code
IF (Tpetra_Have_Kokkos_Refactor)
  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_refactor)
  APPEND_GLOB(HEADERS ${DIR}/*.hpp)
  APPEND_GLOB(SOURCES ${DIR}/*.cpp)
  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})
ENDIF()

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_config.h )
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_ETIHelperMacros.h)

#
# C) Define the targets for package's library(s)
#

# CUDA files
IF (KokkosClassic_ENABLE_Thrust OR
    (KokkosClassic_ENABLE_KokkosCompat AND Kokkos_ENABLE_Cuda))
  FILE(GLOB CUDA_SOURCES *.cu)
ENDIF()

# With the new Kokkos, it is impossible to separate CUDA code from regular
# source code to be compiled with NVCC.  Therefore we compile everything
# with NVCC by including each cpp file in a cmake-generated .cu file.
# We only turn this on if KokkosCompat and Kokkos-Cuda are enabled.  We check
# for KokkosCompat instead of KokkosRefactor because the new Kokkos wrapper
# nodes are enabled when KokkosCompat is enabled, regardless of KokkosRefactor.
IF (KokkosClassic_ENABLE_KokkosCompat AND Kokkos_ENABLE_Cuda)
  FOREACH(CPP_FILE ${SOURCES})
    # Strip .cpp and leading path
    STRING(REGEX REPLACE ".*/(.+).cpp" "\\1" BASE_FILE "${CPP_FILE}")

    # Create .cu file
    SET(CU_FILE "${CMAKE_CURRENT_BINARY_DIR}/${BASE_FILE}_CU.cu")
    CONFIGURE_FILE(Tpetra_nvcc_cufile_template.cu.in ${CU_FILE})

    # Add .cu file to list of nvcc sources
    APPEND_SET(CUDA_SOURCES ${CU_FILE})
  ENDFOREACH()

  # clear SOURCES
  SET(SOURCES "")
ENDIF()

# CUDA library
IF (KokkosClassic_ENABLE_Thrust OR
    (KokkosClassic_ENABLE_KokkosCompat AND Kokkos_ENABLE_Cuda))
  # Add MPI include directory to NVCC flags if MPI enabled
    IF( TPL_ENABLE_MPI )
      CUDA_INCLUDE_DIRECTORIES( ${MPI_BASE_DIR}/include )
    ENDIF()
  PACKAGE_ADD_CUDA_LIBRARY(tpetra_cuda SOURCES ${CUDA_SOURCES})
ENDIF()

TRIBITS_ADD_LIBRARY(
  tpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  ADDED_LIB_TARGET_NAME_OUT TPETRA_LIBNAME
  )

# We need to set the linker language explicitly here for CUDA builds.
SET_PROPERTY(
  TARGET ${TPETRA_LIBNAME}
  APPEND PROPERTY LINKER_LANGUAGE CXX
  )
