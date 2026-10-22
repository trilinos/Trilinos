INCLUDE(PythonTest)

PYTHON_ADD_TEST(test_environment NUM_MPI_PROCS 1,4)
PYTHON_ADD_TEST(test_exit_after_finalize NUM_MPI_PROCS 1,4)
PYTHON_ADD_TEST(test_exit_with_live_mesh NUM_MPI_PROCS 1,4)
PYTHON_ADD_TEST(test_gradient_filter)
PYTHON_ADD_TEST(test_log_redirecter)
PYTHON_ADD_TEST(test_refinement)
PYTHON_ADD_TEST(test_smooth_mesh)
PYTHON_ADD_TEST(test_surface)

