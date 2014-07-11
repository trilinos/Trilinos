INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsCopyFilesToBinaryDir)

SET(HEADERS "")
SET(SOURCES "")


SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

# With the new Kokkos, it is impossible to separate CUDA code from regular
# source code to be compiled with NVCC.  Therefore we compile everything
# with NVCC by including each cpp file in a cmake-generated .cu file.
# We only turn this on if KokkosCompat and Kokkos-Cuda are enabled.  We check
# for KokkosCompat instead of KokkosRefactor because the new Kokkos wrapper
# nodes are enabled when KokkosCompat is enabled, regardless of KokkosRefactor.
IF (KokkosClassic_ENABLE_KokkosCompat AND Kokkos_ENABLE_Cuda)
  CUDA_INCLUDE_DIRECTORIES( ${Tpetra_INCLUDE_DIRS} )
  FOREACH(CPP_FILE ${SOURCES})
    # Strip .cpp and leading path
    STRING(REGEX REPLACE ".*/(.+).cpp" "\\1" BASE_FILE "${CPP_FILE}")

    # Create .cu file
    SET(CU_FILE "${CMAKE_CURRENT_BINARY_DIR}/${BASE_FILE}_CU.cu")
    CONFIGURE_FILE(${DIR}/../../src/Tpetra_nvcc_cufile_template.cu.in ${CU_FILE})

    # Add .cu file to list of nvcc sources
    APPEND_SET(CUDA_SOURCES ${CU_FILE})
  ENDFOREACH()

  # clear SOURCES
  SET(SOURCES "")
  CUDA_COMPILE(CUDA_OBJECTS ${CUDA_SOURCES})
ENDIF()

# The performance test currently depends on KokkosCore for its timers.
# Otherwise, it does not use anything else in KokkosCore directly.
IF (Trilinos_ENABLE_KokkosCore)
  TRIBITS_ADD_EXECUTABLE(
    Performance-CGSolve
    SOURCES
      ${SOURCES} ${CUDA_OBJECTS}
    COMM mpi
  )

#  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyXmlForTpetraPerfCgSolve
#    SOURCE_DIR ${Tpetra_SOURCE_DIR}/test/
#    SOURCE_FILES Tpetra_PerformanceTests.xml
#    EXEDEPS Performance-CGSolve
#  )
ENDIF ()


