INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
 Utst_ioel
 SOURCES Utst_ioel.C
 DEPLIBS Ionit
)

PACKAGE_ADD_TEST(
	Utst_ioel
	NAME Utst_ioel
	NUM_MPI_PROCS 1
)

IF (${PACKAGE_NAME}_ENABLE_SEACASExodus AND ${PACKAGE_NAME}_ENABLE_SEACASNemesis)
  PACKAGE_ADD_EXECUTABLE(
   Utst_superelement
   SOURCES Utst_superelement.C
   DEPLIBS Ionit Ioex
  )

  PACKAGE_ADD_TEST(
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
