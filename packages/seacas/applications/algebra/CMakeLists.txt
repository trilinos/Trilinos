INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(Algebra)

TRIBITS_INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.f)
TRIBITS_ADD_EXECUTABLE(
  algebra
  NOEXEPREFIX
  NOEXESUFFIX
  LINKER_LANGUAGE Fortran
  SOURCES ${SOURCES}
  COMM serial mpi
  INSTALLABLE
  )

TRIBITS_SUBPACKAGE_POSTPROCESS()

if (${CMAKE_PROJECT_NAME}_ENABLE_SEACASExodiff)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/cube_nvar.inp cube_nvar.inp COPYONLY)
TRIBITS_ADD_ADVANCED_TEST(algebra_add_nodal_val
  TEST_0
    EXEC run_algebra.py
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/test/cube_with_ns_ss.g cube_nvar.g cube_nvar.inp
    DIRECTORY ../../scripts
    NOEXEPREFIX NOEXESUFFIX
  TEST_1
    EXEC exodiff
    ARGS -show_all_diffs -pedantic -match_by_name ${CMAKE_CURRENT_SOURCE_DIR}/test/cube_nvar.gold.g cube_nvar.g
    DIRECTORY ../exodiff
    NOEXEPREFIX NOEXESUFFIX
  COMM serial
  XHOSTTYPE Windows
)

endif()
