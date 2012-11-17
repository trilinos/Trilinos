INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  EditCrsMatrix
  SOURCES EditCrsMatrix.cpp
  DEPLIBS xpetra
  )

IF (${PACKAGE_NAME}_ENABLE_Epetra)

	TRIBITS_ADD_TEST(
	  EditCrsMatrix
          NAME "EditCrsMatrixEpetra"
          ARGS "--linAlgebra=Epetra"
          COMM serial mpi
	  )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra)

	TRIBITS_ADD_TEST(
	  EditCrsMatrix
          NAME "EditCrsMatrixTpetra"
          ARGS "--linAlgebra=Tpetra"
          COMM serial mpi
	  )

ENDIF()
