TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../test/unit_tests)

IF (${PACKAGE_NAME}_ENABLE_Belos AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_ADD_EXECUTABLE(
    HierarchicalDriver
    SOURCES main.cpp
    COMM mpi
    )

  file(GLOB FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 1d-*/*)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(caglusa_cp
    SOURCE_FILES
    belos.xml muelu.xml auxiliary.xml ${FILES}
    )

  TRIBITS_ADD_TEST(
    HierarchicalDriver
    NAME "HierarchicalDriver-fractionalLaplacian1D"
    ARGS
    "--xmlHierarchical=1d-mm/hierarchical.xml --xmlProblem=1d-mm/problem.xml"
    "--xmlHierarchical=1d-binary/hierarchical.xml --xmlProblem=1d-binary/problem.xml"
    COMM mpi
    NUM_MPI_PROCS 4
    )

ENDIF()
