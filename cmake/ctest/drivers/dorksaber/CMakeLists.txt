TRILINOS_DRIVER_SETUP()

TRILINOS_DRIVER_ADD_DASHBOARD(
  SERIAL_RELEASE_MueLu_Matlab
  ctest_linux_nightly_serial_release_muelu_matlab_dorksaber.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)


TRILINOS_DRIVER_ADD_DASHBOARD(
  MPI_RELEASE_MueLu_KokkosKernels_Experimental
  ctest_linux_nightly_mpi_release_kokkoskernels_experimental_dorksaber.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 330
)

# Refactor (kokkos) tests
TRILINOS_DRIVER_ADD_DASHBOARD(
  OPENMPI_1.10.0_DEBUG_MueLu_Kokkos_Refactor_OpenMP
  ctest_linux_nightly_mpi_debug_muelu_kokkos_refactor_openmp_dorksaber.cmake
  CTEST_INSTALLER_TYPE release
  RUN_SERIAL
  TIMEOUT_MINUTES 420
  )



TRILINOS_ADD_REQUIRED_CMAKE_INSTALLS()
