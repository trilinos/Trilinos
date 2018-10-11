TRIBITS_ADD_EXECUTABLE(
    stratimikos_frosch
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(StratimikosCopyFiles
    DEST_FILES stratimikos_ParameterList.xml stratimikos_ParameterList.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS  stratimikos_frosch
)

#FR 17/07/2018: Adding Tests
TRIBITS_ADD_TEST(
    stratimikos_frosch
    NAME test_stratimikos_frosch
    ARGS "--DIM=2"
    COMM mpi
    NUM_MPI_PROCS 4
)


