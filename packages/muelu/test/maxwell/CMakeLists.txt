# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Tpetra
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

IF (${PACKAGE_NAME}_ENABLE_Tpetra)
  IF (${PACKAGE_NAME}_ENABLE_Belos)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Maxwell3D
      SOURCES Maxwell3D.cpp
      COMM serial
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Maxwell3D
      SOURCES Maxwell3D.cpp
      NUM_MPI_PROCS 4
      COMM mpi
    )

    IF(HAVE_MUELU_COMPLEX)

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        MaxwellComplex3D
        SOURCES MaxwellComplex3D.cpp
        COMM serial
      )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        MaxwellComplex3D
        SOURCES MaxwellComplex3D.cpp
        NUM_MPI_PROCS 4
        COMM mpi
      )

    ENDIF()

    TRIBITS_COPY_FILES_TO_BINARY_DIR(Maxwell_cp
    SOURCE_FILES M0.txt M1.txt S.txt D0.txt coords.txt
    )
  ENDIF()
ENDIF()

