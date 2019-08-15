# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Belos
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

IF (${PACKAGE_NAME}_ENABLE_Belos AND ${PACKAGE_NAME}_ENABLE_Amesos2 AND ${PACKAGE_NAME}_ENABLE_Galeri)

  TRIBITS_ADD_EXECUTABLE(
    Structured
    SOURCES Driver_Structured.cpp
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_ENABLE_Tpetra)
    TRIBITS_ADD_TEST(
      Structured
      NAME "Structured_Laplace2D_Tpetra"
      ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Laplace2D"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    TRIBITS_ADD_TEST(
      Structured
      NAME "Structured_Laplace2D_Shift_Tpetra"
      ARGS "--linAlgebra=Tpetra --xml=structured_1dof_shift.xml --matrixType=Laplace2D"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    TRIBITS_ADD_TEST(
      Structured
      NAME "Structured_Elasticity3D_Tpetra"
      ARGS "--linAlgebra=Tpetra --xml=structured_3dof.xml --matrixType=Elasticity3D --nx=25 --ny=25 --nz=25"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    TRIBITS_ADD_TEST(
      Structured
      NAME "Structured_Line_Tpetra"
      ARGS "--linAlgebra=Tpetra --xml=structured_interp_line.xml --matrixType=Laplace3D --nx=25 --ny=25 --nz=25"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    IF (${PACKAGE_NAME}_ENABLE_Kokkos_Refactor)
      TRIBITS_ADD_TEST(
        Structured
        NAME "Structured_Interp_Laplace2D_kokkos"
        ARGS "--linAlgebra=Tpetra --xml=structured_interp_kokkos.xml --matrixType=Laplace2D --nx=25 --ny=25"
        COMM serial mpi
        NUM_MPI_PROCS 4
        )

      TRIBITS_ADD_TEST(
        Structured
        NAME "Structured_Interp_SA_Laplace2D_kokkos"
        ARGS "--linAlgebra=Tpetra --xml=structured_interp_sa_kokkos.xml --matrixType=Laplace2D --nx=25 --ny=25"
        COMM serial mpi
        NUM_MPI_PROCS 4
        )

      TRIBITS_ADD_TEST(
        Structured
        NAME "Structured_Tentative_Laplace2D_kokkos"
        ARGS "--linAlgebra=Tpetra --xml=structured_tentative_kokkos.xml --matrixType=Laplace2D --nx=25 --ny=25"
        COMM serial mpi
        NUM_MPI_PROCS 4
        )
    ENDIF()

  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Epetra AND (NOT ${PACKAGE_NAME}_ENABLE_Kokkos_Refactor OR NOT ${PACKAGE_NAME}_ENABLE_Kokkos_Refactor_Use_By_Default))
    TRIBITS_ADD_TEST(
      Structured
      NAME "Structured_Laplace2D_Epetra"
      ARGS "--linAlgebra=Epetra --xml=structured_1dof.xml --matrixType=Laplace2D"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    TRIBITS_ADD_TEST(
      Structured
      NAME "Structured_Laplace2D_Shift_Epetra"
      ARGS "--linAlgebra=Epetra --xml=structured_1dof_shift.xml --matrixType=Laplace2D"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    TRIBITS_ADD_TEST(
      Structured
      NAME "Structured_Elasticity3D_Epetra"
      ARGS "--linAlgebra=Epetra --xml=structured_3dof.xml --matrixType=Elasticity3D --nx=25 --ny=25 --nz=25"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    TRIBITS_ADD_TEST(
      Structured
      NAME "Structured_Line_Epetra"
      ARGS "--linAlgebra=Epetra --xml=structured_interp_line.xml --matrixType=Laplace3D --nx=25 --ny=25 --nz=25"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )
  ENDIF()


  TRIBITS_COPY_FILES_TO_BINARY_DIR(Structured_cp
    SOURCE_FILES structured_1dof.xml structured_1dof_shift.xml structured_2dof.xml structured_3dof.xml structured_interp_line.xml structured_interp_kokkos.xml structured_interp_sa_kokkos.xml structured_tentative_kokkos.xml
    )


ENDIF()
