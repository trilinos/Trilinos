
TRIBITS_ADD_EXECUTABLE(
 Utst_ioel
 SOURCES Utst_ioel.C
)

TRIBITS_ADD_TEST(
	Utst_ioel
	NAME Utst_ioel
	NUM_MPI_PROCS 1
)

TRIBITS_ADD_EXECUTABLE(
 Utst_sort
 SOURCES Utst_sort.C
)

TRIBITS_ADD_TEST(
	Utst_sort
	NAME Utst_sort
	NUM_MPI_PROCS 1
)

TRIBITS_ADD_EXECUTABLE(
 Utst_map
 SOURCES Utst_map.C
)

TRIBITS_ADD_TEST(
	Utst_map
	NAME Utst_map
	NUM_MPI_PROCS 1
)

IF (TPL_ENABLE_CGNS)
TRIBITS_ADD_EXECUTABLE(
 Utst_structured_decomp
 SOURCES Utst_structured_decomp.C
)

TRIBITS_ADD_TEST(
	Utst_structured_decomp
	NAME Utst_structured_decomp
	NUM_MPI_PROCS 1
)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_SEACASExodus)
TRIBITS_ADD_EXECUTABLE(
 Utst_superelement
 SOURCES Utst_superelement.C
)

TRIBITS_ADD_TEST(
	Utst_superelement
	NAME Utst_superelement
	NUM_MPI_PROCS 1
	ARGS ${CMAKE_CURRENT_SOURCE_DIR}/cbr2.ncf
)
ENDIF()

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/../"
  "${CMAKE_CURRENT_BINARY_DIR}/../"
)
