# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Belos)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

IF(${PACKAGE_NAME}_ENABLE_Belos)

  TRIBITS_ADD_EXECUTABLE(
    Structured
    SOURCES Driver_Structured.cpp
    COMM serial mpi
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Structured
    NAME "Structured_Laplace2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Laplace2D"
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Structured
    NAME "Structured_Laplace2D_Shift_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof_shift.xml --matrixType=Laplace2D"
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Structured
    NAME "Structured_Elasticity3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_3dof.xml --matrixType=Elasticity3D --nx=25 --ny=25 --nz=25"
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Structured
    NAME "Structured_Line_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_interp_line.xml --matrixType=Laplace3D --nx=25 --ny=25 --nz=25"
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Structured
    NAME "Structured_Laplace3D_SCP_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_scp_1dof.xml --matrixType=Laplace3D --nx=10 --ny=10 --nz=10"
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Structured
    NAME "Structured_Laplace3D_SPARC_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_sparc_1dof.xml --matrixType=Laplace3D --nx=25 --ny=25 --nz=10"
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Structured
    NAME "Structured_Interp_Laplace2D_kokkos"
    ARGS "--linAlgebra=Tpetra --xml=structured_interp_kokkos.xml --matrixType=Laplace2D --nx=25 --ny=25"
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Structured
    NAME "Structured_Interp_SA_Laplace2D_kokkos"
    ARGS "--linAlgebra=Tpetra --xml=structured_interp_sa_kokkos.xml --matrixType=Laplace2D --nx=25 --ny=25"
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Structured
    NAME "Structured_Tentative_Laplace2D_kokkos"
    ARGS "--linAlgebra=Tpetra --xml=structured_tentative_kokkos.xml --matrixType=Laplace2D --nx=25 --ny=25"
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    Structured
    NAME "StructuredRAPFactory_Constant_Laplace2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structuredrapfactory_laplace_const.xml --matrixType=Laplace2D --nx=10 --ny=10 --mx=1 --my=1 --op=setup --check-rap-against-symbolic"
    COMM serial
    PASS_REGULAR_EXPRESSION "RAP matrix comparison passed: global nnz = 64, graph comparison = exact"
  )

  TRIBITS_ADD_TEST(
    Structured
    NAME "StructuredRAPFactory_Constant_Laplace2D_XYSplit_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structuredrapfactory_laplace_const.xml --matrixType=Laplace2D --nx=10 --ny=10 --mx=2 --my=2 --op=setup --check-rap-against-symbolic"
    COMM mpi
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "RAP matrix comparison passed: global nnz = 156, graph comparison = exact"
  )

  TRIBITS_ADD_TEST(
    Structured
    NAME "StructuredRAPFactory_Constant_Laplace3D_ZSplit_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structuredrapfactory_laplace_const.xml --matrixType=Laplace3D --nx=10 --ny=10 --nz=20 --mx=1 --my=1 --mz=2 --op=setup --check-rap-against-symbolic"
    COMM mpi
    NUM_MPI_PROCS 2
    PASS_REGULAR_EXPRESSION "RAP matrix comparison passed: global nnz = 736, graph comparison = exact"
  )

  TRIBITS_ADD_TEST(
    Structured
    NAME "StructuredRAPFactory_Linear_Laplace3D_ZSplit_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structuredrapfactory_laplace_lin.xml --matrixType=Laplace3D --nx=10 --ny=10 --nz=20 --mx=1 --my=1 --mz=2 --op=setup --check-rap-against-symbolic --rap-graph-comparison=prebuilt-superset"
    COMM mpi
    NUM_MPI_PROCS 2
    PASS_REGULAR_EXPRESSION "RAP matrix comparison passed: global nnz = 8704, graph comparison = prebuilt superset, symbolic global nnz = 8232, extra prebuilt entries = 472"
  )

  TRIBITS_ADD_TEST(
    Structured
    NAME "StructuredRAPFactory_Constant_Elasticity2D_XYSplit_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structuredrapfactory_elasticity2d_const.xml --matrixType=Elasticity2D --nx=20 --ny=20 --mx=2 --my=2 --op=setup --check-rap-against-symbolic"
    COMM mpi
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "RAP matrix comparison passed: global nnz = 1936, graph comparison = exact"
  )

  TRIBITS_ADD_TEST(
    Structured
    NAME "StructuredRAPFactory_Constant_Elasticity3D_XYSplit_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structuredrapfactory_elasticity3d_const.xml --matrixType=Elasticity3D --nx=20 --ny=20 --nz=10 --mx=2 --my=2 --mz=1 --op=setup --check-rap-against-symbolic"
    COMM mpi
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "RAP matrix comparison passed: global nnz = 43560, graph comparison = exact"
  )

  TRIBITS_ADD_TEST(
    Structured
    NAME "StructuredRAPFactory_Linear_Elasticity3D_XYSplit_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structuredrapfactory_elasticity3d_lin.xml --matrixType=Elasticity3D --nx=16 --ny=16 --nz=16 --mx=2 --my=2 --mz=1 --op=setup --check-rap-against-symbolic"
    COMM mpi
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "RAP matrix comparison passed: global nnz = 29952, graph comparison = exact"
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    Structured_cp
    SOURCE_FILES
      structured_1dof.xml
      structured_1dof_shift.xml
      structured_2dof.xml
      structured_3dof.xml
      structured_interp_kokkos.xml
      structured_interp_line.xml
      structured_interp_sa_kokkos.xml
      structured_scp_1dof.xml
      structured_sparc_1dof.xml
      structured_tentative_kokkos.xml
      structuredrapfactory_elasticity2d_const.xml
      structuredrapfactory_elasticity3d_const.xml
      structuredrapfactory_elasticity3d_lin.xml
      structuredrapfactory_laplace_const.xml
      structuredrapfactory_laplace_lin.xml
  )

ENDIF()
