INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (${PACKAGE_NAME}_ENABLE_Tpetra)  
	TRIBITS_ADD_EXECUTABLE_AND_TEST(
	  Laplace1D-Xpetra
	  SOURCES Laplace1D-Xpetra.cpp
	  DEPLIBS xpetra
	  )


	TRIBITS_ADD_EXECUTABLE_AND_TEST(
	  Laplace1D-Tpetra
	  SOURCES Laplace1D-Tpetra.cpp
	  DEPLIBS xpetra
	  )
ENDIF()
