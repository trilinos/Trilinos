# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Tpetra
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_Galeri
  )

INCLUDE(TribitsAddExecutable)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

IF (${PACKAGE_NAME}_ENABLE_Galeri)

  IF(HAVE_MUELU_COMPLEX)

    TRIBITS_ADD_EXECUTABLE(
      Helmholtz1D
      SOURCES Helmholtz1D.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra muelu-adapters
    )

    TRIBITS_ADD_EXECUTABLE(
      Helmholtz2D
      SOURCES Helmholtz2D.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra muelu-adapters
    )

    TRIBITS_ADD_EXECUTABLE(
      Helmholtz3D
      SOURCES Helmholtz3D.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra muelu-adapters
    )

    TRIBITS_ADD_EXECUTABLE(
      HelmholtzFEM2D
      SOURCES HelmholtzFEM2D.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra muelu-adapters
    )

    TRIBITS_ADD_EXECUTABLE(
      HelmholtzFEM3D
      SOURCES HelmholtzFEM3D.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra muelu-adapters
    )

  ENDIF()

ENDIF() # Galeri

