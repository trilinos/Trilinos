
TRIBITS_ADD_EXECUTABLE(
 Utst_textmesh
 SOURCES unitMain.C UnitTestTextMesh.C
)

TRIBITS_ADD_TEST(
        Utst_textmesh
        NAME Utst_textmesh
        NUM_MPI_PROCS 1
)

TRIBITS_ADD_TEST(
        Utst_textmesh
        NAME Utst_textmesh
        NUM_MPI_PROCS 2
)

TRIBITS_ADD_TEST(
        Utst_textmesh
        NAME Utst_textmesh
        NUM_MPI_PROCS 3
)

TRIBITS_ADD_TEST(
        Utst_textmesh
        NAME Utst_textmesh
        NUM_MPI_PROCS 4
)

TRIBITS_INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/../"
  "${CMAKE_CURRENT_BINARY_DIR}/../"
)
