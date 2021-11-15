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
    "--xml=hierarchical-1d-mm-global.xml"
    "--xml=hierarchical-1d-mm-local.xml"
    "--xml=hierarchical-1d-binary-global.xml"
    "--xml=hierarchical-1d-binary-local.xml"
    COMM mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(caglusa_cp
    SOURCE_FILES
    hierarchical-1d-mm-global.xml hierarchical-1d-mm-local.xml hierarchical-1d-binary-global.xml hierarchical-1d-binary-local.xml muelu.xml aux.xml
    )

  ADD_SUBDIRECTORY(1d-mm-global)
  ADD_SUBDIRECTORY(1d-mm-local)
  ADD_SUBDIRECTORY(1d-binary-global)
  ADD_SUBDIRECTORY(1d-binary-local)

ENDIF()
