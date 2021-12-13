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
    COMM mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(caglusa_cp
    SOURCE_FILES
    hierarchical.xml muelu.xml aux.xml
    rowmap.mtx coeffmap.mtx ghosted_coeffmap.mtx
    Anear.mtx basismatrix.mtx kernelApproximations.mtx
    transfer0.mtx transfer1.mtx transfer2.mtx transfer3.mtx transfer4.mtx
    X_ex.mtx RHS.mtx
    aux.mtx coords.mtx
    )

ENDIF()
