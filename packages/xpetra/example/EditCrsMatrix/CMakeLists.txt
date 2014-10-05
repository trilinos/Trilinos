INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (NOT Xpetra_Epetra_NO_32BIT_GLOBAL_INDICES)
  TRIBITS_ADD_EXECUTABLE(
    EditCrsMatrix
    SOURCES EditCrsMatrix.cpp
    DEPLIBS xpetra
    )
ENDIF()


IF (NOT Xpetra_Epetra_NO_64BIT_GLOBAL_INDICES)
  TRIBITS_ADD_EXECUTABLE(
    EditCrsMatrix_LL
    SOURCES EditCrsMatrix.cpp
    DEPLIBS xpetra
    DEFINES -DXPETRA_TEST_USE_LONGLONG_GO
    )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Epetra)

  IF (NOT Xpetra_Epetra_NO_32BIT_GLOBAL_INDICES)
	TRIBITS_ADD_TEST(
	  EditCrsMatrix
          NAME "EditCrsMatrixEpetra"
          ARGS "--linAlgebra=Epetra"
          COMM serial mpi
	  )
  ENDIF()
  IF (NOT Xpetra_Epetra_NO_64BIT_GLOBAL_INDICES)
	TRIBITS_ADD_TEST(
      EditCrsMatrix_LL
          NAME "EditCrsMatrixEpetra"
          ARGS "--linAlgebra=Epetra"
          COMM serial mpi
	  )
  ENDIF()

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra)

  IF (NOT Xpetra_Epetra_NO_32BIT_GLOBAL_INDICES)
	TRIBITS_ADD_TEST(
	  EditCrsMatrix
          NAME "EditCrsMatrixTpetra"
          ARGS "--linAlgebra=Tpetra"
          COMM serial mpi
	  )
  ENDIF()
  IF (NOT Xpetra_Epetra_NO_64BIT_GLOBAL_INDICES)
	TRIBITS_ADD_TEST(
      EditCrsMatrix_LL
          NAME "EditCrsMatrixTpetra"
          ARGS "--linAlgebra=Tpetra"
          COMM serial mpi
	  )
  ENDIF()

ENDIF()
