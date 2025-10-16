# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos2 ${PACKAGE_NAME}_ENABLE_Belos)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

IF(${PACKAGE_NAME}_ENABLE_Belos AND ${PACKAGE_NAME}_ENABLE_Amesos2)

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

  IF(${PACKAGE_NAME}_ENABLE_Epetra)
    MUELU_ADD_SERIAL_AND_MPI_TEST(
      Structured
      NAME "Structured_Laplace2D_Epetra"
      ARGS "--linAlgebra=Epetra --xml=structured_1dof.xml --matrixType=Laplace2D"
      COMM serial mpi
      NUM_MPI_PROCS 4
    )

    MUELU_ADD_SERIAL_AND_MPI_TEST(
      Structured
      NAME "Structured_Laplace2D_Shift_Epetra"
      ARGS "--linAlgebra=Epetra --xml=structured_1dof_shift.xml --matrixType=Laplace2D"
      COMM serial mpi
      NUM_MPI_PROCS 4
    )

    MUELU_ADD_SERIAL_AND_MPI_TEST(
      Structured
      NAME "Structured_Elasticity3D_Epetra"
      ARGS "--linAlgebra=Epetra --xml=structured_3dof.xml --matrixType=Elasticity3D --nx=25 --ny=25 --nz=25"
      COMM serial mpi
      NUM_MPI_PROCS 4
    )

    MUELU_ADD_SERIAL_AND_MPI_TEST(
      Structured
      NAME "Structured_Line_Epetra"
      ARGS "--linAlgebra=Epetra --xml=structured_interp_line.xml --matrixType=Laplace3D --nx=25 --ny=25 --nz=25"
      COMM serial mpi
      NUM_MPI_PROCS 4
    )
  ENDIF()

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
  )

ENDIF()
