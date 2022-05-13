INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../test/unit_tests)

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Belos AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_ADD_EXECUTABLE(
    HierarchicalDriver
    SOURCES main.cpp
    COMM mpi
    )

  file(GLOB FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 1d-*/*)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(caglusa_cp
    SOURCE_FILES
    muelu.xml aux.xml ${FILES}
    )

  TRIBITS_ADD_TEST(
    HierarchicalDriver
    NAME "HierarchicalDriver-fractionalLaplacian1D"
    ARGS
    "--xml=1d-mm-global/hierarchical.xml"
    "--xml=1d-mm-local/hierarchical.xml"
    "--xml=1d-binary-global/hierarchical.xml"
    "--xml=1d-binary-local/hierarchical.xml"
    COMM mpi
    NUM_MPI_PROCS 4
    )

ENDIF()
