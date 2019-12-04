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

IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
localpartitionofunitybasis
NAME test_localpartitionofunitybasis_EPETRA
COMM mpi
ARGS "--USEEPETRA"
NUM_MPI_PROCS 1
)
ENDIF()

TRIBITS_ADD_TEST(
localpartitionofunitybasis
NAME test_localpartitionofunitybasis_TPETRA
ARGS "--USETPETRA"
COMM mpi
NUM_MPI_PROCS 1
)
