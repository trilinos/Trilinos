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
    SOURCE_FILES amg_1dof.xml
      structured_1dof.xml
      structured_1dof_3level.xml
      structured_1dof-complex.xml
      structured_linear_1dof.xml
      structured_unstructured_1dof.xml
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_1.log"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_4.log"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_AMG_CoarseSolver_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof_3level.xml --matrixType=Star2D --nx=50 --ny=50 --smootherIts=2 --amg-coarse-solver --coarseAmgXml=amg_1dof.xml"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Star2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_linear_1.log"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Star2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_linear_4.log"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Brick3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Brick3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Brick3D_1.log"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Brick3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Brick3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Brick3D_4.log"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Brick3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_1dof.xml --matrixType=Brick3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Brick3D_linear_1.log"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Brick3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_1dof.xml --matrixType=Brick3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Brick2D_linear_4.log"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Unstructured_Laplace3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_unstructured_1dof.xml --matrixType=Laplace3D --nx=9 --ny=9 --nz=5 --smootherIts=2 --unstructured={3}"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

ENDIF()

