
ENABLE_TESTING()

INCLUDE(PyTrilinos2MakeTest)

PyTrilinos2_MAKE_MPI_TEST(CG)
PyTrilinos2_MAKE_MPI_TEST(parameterList)
