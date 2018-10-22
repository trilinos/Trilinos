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
    Maxwell3D
    SOURCES Maxwell3D.cpp
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_ENABLE_Tpetra)
    TRIBITS_ADD_TEST(
      Maxwell3D
      NAME "Maxwell3D-Tpetra"
      ARGS "--linAlgebra=Tpetra"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    IF (${PACKAGE_NAME}_INST_DOUBLE_INT_INT)
      TRIBITS_ADD_TEST(
        Maxwell3D
        NAME "Maxwell3D-Tpetra-Stratimikos"
        ARGS "--linAlgebra=Tpetra --solverName=Stratimikos --instantiation=DOUBLE_INT_INT"
        COMM serial mpi
        NUM_MPI_PROCS 4
        )
    ENDIF()
  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Epetra)
    TRIBITS_ADD_TEST(
      Maxwell3D
      NAME "Maxwell3D-Epetra"
      ARGS "--linAlgebra=Epetra"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )
  ENDIF()


  TRIBITS_COPY_FILES_TO_BINARY_DIR(Maxwell_cp
    SOURCE_FILES M0.mat M1.mat S.mat D0.mat coords.mat Maxwell.xml
    )

  IF(HAVE_MUELU_COMPLEX)
    TRIBITS_COPY_FILES_TO_BINARY_DIR(Maxwell_cp_complex
      SOURCE_FILES M0_complex.mat M1_complex.mat S_complex.mat D0_complex.mat Maxwell_complex.xml
      )
  ENDIF()


ENDIF()
