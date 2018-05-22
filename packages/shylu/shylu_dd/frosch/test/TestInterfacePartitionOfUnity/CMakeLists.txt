TRIBITS_ADD_EXECUTABLE(
	frosch_interfacepartitionofunity
	SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLU_DDFROSchIPOUCopyFiles
    DEST_FILES ParametersIPOU.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS frosch_interfacepartitionofunity
)

TRIBITS_ADD_TEST(
	frosch_interfacepartitionofunity
	NAME test_frosch_frosch_interfacepartitionofunity
	ARGS "--DIM=2 --M=4"
	COMM mpi
	NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
    frosch_interfacepartitionofunity
    NAME test_frosch_interfacepartitionofunity
    ARGS "--DIM=3 --M=4"
    COMM mpi
    NUM_MPI_PROCS 8
)
