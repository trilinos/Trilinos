TRIBITS_ADD_EXECUTABLE(
	test_import
	SOURCES main.cpp
)

IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
TRIBITS_ADD_TEST(
test_import
NAME test_import_EPETRA
NUM_MPI_PROCS 8
ARGS "--USEEPETRA"
)
ENDIF()

TRIBITS_ADD_TEST(
test_import
NAME test_import_TPETRA
NUM_MPI_PROCS 8
ARGS "--USETPETRA"
)
