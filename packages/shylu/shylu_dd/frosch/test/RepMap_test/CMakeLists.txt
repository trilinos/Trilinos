TRIBITS_ADD_EXECUTABLE(
    repmap
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(RepMapCopyFiles
    DEST_FILES xpetra_ParameterList.xml
    DEST_FILES GDSW.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS repmap
)

TRIBITS_ADD_TEST(
    repmap
    NAME test_repmap
    COMM mpi
    NUM_MPI_PROCS 4
)