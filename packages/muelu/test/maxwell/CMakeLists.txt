# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos2
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

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Belos AND ${PACKAGE_NAME}_ENABLE_Amesos2)

# TAW: 10/15/2015  Disabling until bug#6358 is fixed for MueLu tests.
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

# JJH: 8/21/2015  Disabling until LinearSolverFactory is ready.  See bug#6395.
#    IF(HAVE_MUELU_COMPLEX)
#
#      TRIBITS_ADD_EXECUTABLE_AND_TEST(
#        MaxwellComplex3D
#        SOURCES MaxwellComplex3D.cpp
#        COMM serial
#      )
#
#      TRIBITS_ADD_EXECUTABLE_AND_TEST(
#        MaxwellComplex3D
#        SOURCES MaxwellComplex3D.cpp
#        NUM_MPI_PROCS 4
#        COMM mpi
#      )
#
#    ENDIF()

    TRIBITS_COPY_FILES_TO_BINARY_DIR(Maxwell_cp
    SOURCE_FILES M0.txt M1.txt S.txt D0.txt coords.txt
    )
ENDIF()

