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

IF (${PACKAGE_NAME}_ENABLE_Belos AND ${PACKAGE_NAME}_ENABLE_Amesos2)

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
      NAME "Structured_Elasticity3D_Tpetra"
      ARGS "--linAlgebra=Tpetra --xml=structured_3dof.xml --matrixType=Elasticity3D --nx=25 --ny=25 --nz=25"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )
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
      NAME "Structured_Elasticity3D_Epetra"
      ARGS "--linAlgebra=Epetra --xml=structured_3dof.xml --matrixType=Elasticity3D --nx=25 --ny=25 --nz=25"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )
  ENDIF()


  TRIBITS_COPY_FILES_TO_BINARY_DIR(Structured_cp
    SOURCE_FILES structured_1dof.xml structured_2dof.xml structured_3dof.xml
    )


ENDIF()
