INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(Kokkos_ENABLE_SERIAL)
  ADD_SUBDIRECTORY(Serial)
ENDIF()

IF(Kokkos_ENABLE_OPENMP)
  ADD_SUBDIRECTORY(OpenMP)
ENDIF()

IF(Kokkos_ENABLE_CUDA AND Kokkos_ENABLE_CUDA_UVM)
  ADD_SUBDIRECTORY(Cuda)
ENDIF()


TRIBITS_COPY_FILES_TO_BINARY_DIR(IntegrationTestDataCopy
  SOURCE_FILES
    EDGE_integrals.dat
    HEX_integrals.dat
    PYR_integrals.dat
    QUAD_integrals.dat
    TET_integrals.dat
    TRI_integrals.dat
    TRIPRISM_integrals.dat
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/data"
#  EXEDEPS Test_01 #Test_02
)

