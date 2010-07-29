INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
 Utst_ioel
 SOURCES Utst_ioel.C
 DEPLIBS Ionit
 COMM serial
)

PACKAGE_ADD_TEST(
	Utst_ioel
	NAME Utst_ioel
	NUM_MPI_PROCS 1
	COMM serial
)

PACKAGE_ADD_EXECUTABLE(
 Utst_superelement
 SOURCES Utst_superelement.C
 DEPLIBS Ionit
 COMM serial
)

PACKAGE_ADD_TEST(
	Utst_superelement
	NAME Utst_superelement
	NUM_MPI_PROCS 1
	COMM serial
	ARGS ${CMAKE_CURRENT_SOURCE_DIR}/cbr2.ncf
)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/../"
  "${CMAKE_CURRENT_BINARY_DIR}/../"
)
