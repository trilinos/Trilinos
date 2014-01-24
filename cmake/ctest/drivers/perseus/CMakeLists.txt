TRILINOS_DRIVER_SETUP()

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_Tpetra_KokkosRefactor_CudaOpenMP_GCC
  ctest_linux_nightly_mpi_release_tpetra_kokkosrefactor_perseus.cmake 
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 60
)


TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
