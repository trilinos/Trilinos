TRIBITS_ADD_EXECUTABLE(
	localpartitionofunitybasis
	SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(LocalPOUCopyFiles
    DEST_FILES Parameters.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS localpartitionofunitybasis
)

TRIBITS_ADD_TEST(
	localpartitionofunitybasis
	NAME test_localpartitionofunitybasis
    COMM mpi
    NUM_MPI_PROCS 1
)
