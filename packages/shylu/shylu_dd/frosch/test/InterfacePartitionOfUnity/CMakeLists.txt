TRIBITS_ADD_EXECUTABLE(
	interfacepartitionofunity
	SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(IPOUCopyFiles
    DEST_FILES ParametersIPOU.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS interfacepartitionofunity
)

TRIBITS_ADD_TEST(
	interfacepartitionofunity
	NAME test_interfacepartitionofunity
	ARGS "--DIM=2 --M=4"
	COMM mpi
	NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
    interfacepartitionofunity
    NAME test_interfacepartitionofunity
    ARGS "--DIM=3 --M=4"
    COMM mpi
    NUM_MPI_PROCS 8
)
