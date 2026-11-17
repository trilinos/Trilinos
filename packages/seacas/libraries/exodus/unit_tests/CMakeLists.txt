TRIBITS_ADD_EXECUTABLE(
 Utst_exodus_names
 SOURCES UnitTestNameLength.C
)

TRIBITS_ADD_TEST(
        Utst_exodus_names
        NAME Utst_exodus_names
        NUM_MPI_PROCS 1
)
