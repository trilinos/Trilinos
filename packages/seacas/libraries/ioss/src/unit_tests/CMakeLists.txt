
TRIBITS_ADD_EXECUTABLE(
 Utst_textmesh
 SOURCES unitMain.C UnitTestTextMesh.C
)

TRIBITS_ADD_EXECUTABLE(
 Utst_iotm
 SOURCES unitMain.C UnitTestIotmDatabaseIO.C
)

TRIBITS_ADD_EXECUTABLE(
 Utst_blockbatchread
 SOURCES unitMain.C UnitTestElementBlockBatchRead.C
)

TRIBITS_ADD_EXECUTABLE(
 Utst_dynamictopology
 SOURCES unitMain.C UnitTestDynamicTopology.C
)

TRIBITS_ADD_EXECUTABLE(
 Utst_generatedmesh
 SOURCES unitMain.C UnitTestGeneratedMesh.C
)

TRIBITS_ADD_EXECUTABLE(
 Utst_flushhandler
 SOURCES unitMain.C UnitTestFlushHandler.C
)

TRIBITS_ADD_EXECUTABLE(
 Utst_serialparallelfilename
 SOURCES unitMain.C UnitTestSerialParallelFilename.C
)

TRIBITS_ADD_TEST(
        Utst_generatedmesh
        NAME Utst_generatedmesh
        NUM_MPI_PROCS 1
)

TRIBITS_ADD_TEST(
        Utst_dynamictopology
        NAME Utst_dynamictopology
        NUM_MPI_PROCS 1
)

TRIBITS_ADD_TEST(
        Utst_blockbatchread
        NAME Utst_blockbatchread
        NUM_MPI_PROCS 2
        ADDED_TESTS_NAMES_OUT BlockBatchRead2AddedName
)

TRIBITS_ADD_TEST(
        Utst_blockbatchread
        NAME Utst_blockbatchread
        NUM_MPI_PROCS 4
        ADDED_TESTS_NAMES_OUT BlockBatchRead4AddedName
)

TRIBITS_ADD_TEST(
        Utst_blockbatchread
        NAME Utst_blockbatchread
        NUM_MPI_PROCS 8
        ADDED_TESTS_NAMES_OUT BlockBatchRead8AddedName
)

# The blockbatchread tests write and read a generated mesh file in the
# same directory; serialize them so they cannot interfere with each other.
set(BlockBatchReadAddedNames
    ${BlockBatchRead2AddedName} ${BlockBatchRead4AddedName} ${BlockBatchRead8AddedName})
if (BlockBatchReadAddedNames)
  set_tests_properties(${BlockBatchReadAddedNames}
    PROPERTIES RESOURCE_LOCK Utst_blockbatchread_mesh_file)
endif()

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

TRIBITS_ADD_TEST(
        Utst_iotm
        NAME Utst_iotm
        NUM_MPI_PROCS 1
)

TRIBITS_ADD_TEST(
        Utst_iotm
        NAME Utst_iotm
        NUM_MPI_PROCS 2
)

TRIBITS_ADD_TEST(
        Utst_flushhandler
        NAME  Utst_flushhandler
        NUM_MPI_PROCS 1
)

TRIBITS_ADD_TEST(
        Utst_flushhandler
        NAME  Utst_flushhandler
        NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
        Utst_serialparallelfilename
        NAME Utst_serialparallelfilename
        NUM_MPI_PROCS 1
)

TRIBITS_INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/../"
  "${CMAKE_CURRENT_BINARY_DIR}/../"
)
