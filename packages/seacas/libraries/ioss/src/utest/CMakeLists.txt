INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
 Utst_ioel
 SOURCES Utst_ioel.C
)

TRIBITS_ADD_TEST(
	Utst_ioel
	NAME Utst_ioel
	NUM_MPI_PROCS 1
)

IF (${PACKAGE_NAME}_ENABLE_EXODUS)
TRIBITS_ADD_EXECUTABLE(
 Utst_superelement
 SOURCES Utst_superelement.C
)

TRIBITS_ADD_TEST(
	Utst_superelement
	NAME Utst_superelement
	NUM_MPI_PROCS 1
	COMM serial
	ARGS ${CMAKE_CURRENT_SOURCE_DIR}/cbr2.ncf
)
ENDIF()

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/../"
  "${CMAKE_CURRENT_BINARY_DIR}/../"
)
