INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  EditCrsMatrix.cpp
  SOURCES EditCrsMatrix.cpp
  DEPLIBS xpetra
  )

IF (${PACKAGE_NAME}_ENABLE_Epetra)  

	TRIBITS_ADD_TEST(
	  EditCrsMatrix.cpp
	  SOURCES EditCrsMatrix.cpp
          ARGS "--linAlgebra=0"
          COMM serial mpi
	  )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra)  

	TRIBITS_ADD_TEST(
	  EditCrsMatrix.cpp
	  SOURCES EditCrsMatrix.cpp
          ARGS "--linAlgebra=1"
          COMM serial mpi
	  )

ENDIF()
