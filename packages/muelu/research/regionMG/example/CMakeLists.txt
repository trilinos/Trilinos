TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../test/unit_tests)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src)

# This test requires Tpetra and Amesos2, so it's only included if Tpetra is enabled
IF(${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_ADD_EXECUTABLE(
    RegionExampleDriver
    SOURCES ../test/structured/Driver_Structured_Regions.cpp # Reuse the source code form the test directory
    COMM serial mpi
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    StructuredRegion_example_cp
    SOURCE_FILES
      amg_1dof.xml
      elasticity_3d.xml
      poisson_3d.xml
      structured_unstructured_1dof.xml
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    RegionExampleDriver
    NAME "Structured_Region_Example_Elasticity3D"
    ARGS "--linAlgebra=Tpetra --xml=elasticity_3d.xml --matrixType=Elasticity3D --nx=28 --ny=28 --nz=28 --smootherIts=2"
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    RegionExampleDriver
    NAME "Structured_Region_Example_Poisson3D_27pt_stencil"
    ARGS "--linAlgebra=Tpetra --xml=poisson_3d.xml --matrixType=Brick3D --nx=28 --ny=28 --nz=28--smootherIts=2"
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    RegionExampleDriver
    NAME "Structured_Region_Example_Poisson3D_27pt_stencil_AMG"
    ARGS "--linAlgebra=Tpetra --xml=poisson_3d.xml --matrixType=Brick3D --nx=50 --ny=50 --nz=50 --smootherIts=2 --coarseSolverType=amg --coarseAmgXml=amg_1dof.xml"
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    RegionExampleDriver
    NAME "Structured_Unstructured_Region_Example_Poisson3D_27pt_stencil"
    ARGS "--linAlgebra=Tpetra --xml=structured_unstructured_1dof.xml --matrixType=Brick3D --nx=28 --ny=28 --nz=28 --smootherIts=2 --unstructured={1,2}"
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

ENDIF()
