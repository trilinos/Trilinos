# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_Zoltan2Core
)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

# mfh 18 Aug 2015: All of the Helmholtz examples require complex arithmetic.
IF(${PACKAGE_NAME}_ENABLE_Zoltan2Core AND HAVE_TPETRA_COMPLEX_DOUBLE)

  TRIBITS_ADD_EXECUTABLE(
    Helmholtz1D
    SOURCES Helmholtz1D.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    Helmholtz2D
    SOURCES Helmholtz2D.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    Helmholtz3D
    SOURCES Helmholtz3D.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    HelmholtzFEM2D
    SOURCES HelmholtzFEM2D.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    HelmholtzFEM3D
    SOURCES HelmholtzFEM3D.cpp
    COMM serial mpi
  )

  # TRIBITS_ADD_TEST(
  #   CreateOperator
  #   NAME "CreateOperatorEpetra"
  #   ARGS "--linAlgebra=Epetra"
  #   STANDARD_PASS_OUTPUT
  #   NUM_MPI_PROCS 1
  #   COMM mpi serial
  #   )

  TRIBITS_ADD_TEST(
    Helmholtz1D
    NAME "Helmholtz1DSerial"
    COMM serial
  )

  TRIBITS_ADD_TEST(
    Helmholtz2D
    NAME "Helmholtz2DSerial"
    CATEGORIES NIGHTLY
    COMM serial
  )

  TRIBITS_ADD_TEST(
    Helmholtz3D
    NAME "Helmholtz3DSerial"
    COMM serial
  )

  TRIBITS_ADD_TEST(
    HelmholtzFEM2D
    NAME "HelmholtzFEM2DSerial"
    COMM serial
  )

  TRIBITS_ADD_TEST(
    HelmholtzFEM3D
    NAME "HelmholtzFEM3DSerial"
    CATEGORIES NIGHTLY
    COMM serial
  )

  TRIBITS_ADD_TEST(
    Helmholtz1D
    NAME "Helmholtz1DParallel"
    NUM_MPI_PROCS 4
    COMM mpi
  )

  TRIBITS_ADD_TEST(
    Helmholtz2D
    NAME "Helmholtz2DParallel"
    CATEGORIES CONTINUOUS
    NUM_MPI_PROCS 4
    COMM mpi
  )

  TRIBITS_ADD_TEST(
    Helmholtz3D
    NAME "Helmholtz3DParallel"
    NUM_MPI_PROCS 4
    COMM mpi
  )

  TRIBITS_ADD_TEST(
    HelmholtzFEM2D
    NAME "HelmholtzFEM2DParallel"
    NUM_MPI_PROCS 4
    COMM mpi
  )

  TRIBITS_ADD_TEST(
    HelmholtzFEM3D
    NAME "HelmholtzFEM3DParallel"
    NUM_MPI_PROCS 4
    COMM mpi
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Helmholtz_cp SOURCE_FILES helm2D.inp helm3D.inp)

ENDIF()
