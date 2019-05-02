INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../../test/unit_tests)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../mmayr/composite_to_regions/src)

# This test requires Tpetra, so it's only included if Tpetra is enabled
IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_ADD_EXECUTABLE(
    StructuredRegionDriver
    SOURCES Driver_Structured_Regions.cpp
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(StructuredRegion_cp
    SOURCE_FILES structured_1dof.xml structured_linear_1dof.xml
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Star2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_1dof.xml --matrixType=Star2D --nx=10 --ny=10"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Star2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_1dof.xml --matrixType=Star2D --nx=10 --ny=10"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

ENDIF()

