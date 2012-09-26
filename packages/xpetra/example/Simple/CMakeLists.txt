INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../gallery)

IF (${PACKAGE_NAME}_ENABLE_Tpetra)  
	TRIBITS_ADD_EXECUTABLE_AND_TEST(
	  Xpetra_Laplace1D-Xpetra
	  SOURCES Laplace1D-Xpetra.cpp
	  DEPLIBS xpetra
	  )


	TRIBITS_ADD_EXECUTABLE_AND_TEST(
	  Xpetra_Laplace1D-Tpetra
	  SOURCES Laplace1D-Tpetra.cpp
	  DEPLIBS xpetra
	  )
ENDIF()
