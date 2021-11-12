INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../test/unit_tests)

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Belos AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_ADD_EXECUTABLE(
    HierarchicalDriver
    SOURCES main.cpp
    COMM mpi
    )

  TRIBITS_ADD_TEST(
    HierarchicalDriver
    NAME "HierarchicalDriver-fractionalLaplacian1D"
    ARGS
    "--xml=hierarchical-1d-mm.xml"
    "--xml=hierarchical-1d-binary.xml"
    COMM mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(caglusa_cp
    SOURCE_FILES
    hierarchical-1d-mm.xml hierarchical-1d-binary.xml muelu.xml aux.xml
    )

  ADD_SUBDIRECTORY(1d-mm)
  ADD_SUBDIRECTORY(1d-binary)

ENDIF()
