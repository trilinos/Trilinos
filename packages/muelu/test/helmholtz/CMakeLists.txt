# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Tpetra
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_Galeri
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

IF (${PACKAGE_NAME}_ENABLE_Galeri)

  IF(HAVE_MUELU_COMPLEX)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Helmholtz1D
      SOURCES Helmholtz1D.cpp
      COMM serial
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Helmholtz2D
      SOURCES Helmholtz2D.cpp
      COMM serial
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Helmholtz3D
      SOURCES Helmholtz3D.cpp
      COMM serial
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      HelmholtzFEM2D
      SOURCES HelmholtzFEM2D.cpp
      COMM serial
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      HelmholtzFEM3D
      SOURCES HelmholtzFEM3D.cpp
      COMM serial
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Helmholtz1D
      SOURCES Helmholtz1D.cpp
      NUM_MPI_PROCS 4
      COMM mpi
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Helmholtz2D
      SOURCES Helmholtz2D.cpp
      NUM_MPI_PROCS 4
      COMM mpi
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Helmholtz3D
      SOURCES Helmholtz3D.cpp
      NUM_MPI_PROCS 4
      COMM mpi
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      HelmholtzFEM2D
      SOURCES HelmholtzFEM2D.cpp
      NUM_MPI_PROCS 4
      COMM mpi
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      HelmholtzFEM3D
      SOURCES HelmholtzFEM3D.cpp
      NUM_MPI_PROCS 4
      COMM mpi
    )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(Helmholtz_cp
    SOURCE_FILES helm2D.inp helm3D.inp
    )

  ENDIF()

ENDIF() # Galeri

