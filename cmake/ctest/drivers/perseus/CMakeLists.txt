TRILINOS_DRIVER_SETUP()

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_Tpetra_KokkosRefactor_GCC
  ctest_linux_nightly_mpi_release_tpetra_kokkosrefactor_perseus.cmake 
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 60
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_Tpetra_KokkosRefactor_ETI_CUDA_OPENMP_HWLOC_GCC
  ctest_linux_nightly_mpi_release_tpetra_kokkosrefactor_ETI_cuda_openmp_hwloc_gcc_perseus.cmake 
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 60
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_Tpetra_KokkosRefactor_CUDA_OPENMP_HWLOC_GCC
  ctest_linux_nightly_mpi_release_tpetra_kokkosrefactor_cuda_openmp_hwloc_gcc_perseus.cmake 
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 60
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_Downstream_KokkosRefactor_ETI_OPENMP_HWLOC_GCC
  ctest_linux_nightly_mpi_release_downstream_kokkosrefactor_ETI_openmp_hwloc_gcc_perseus.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 180
)

TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_Downstream_KokkosRefactor_ETI_SERIAL_Intel
  ctest_linux_nightly_mpi_release_downstream_kokkosrefactor_ETI_serial_intel_perseus.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 180
)


TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
