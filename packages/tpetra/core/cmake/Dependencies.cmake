TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_PACKAGES Teuchos Kokkos TeuchosKokkosCompat TeuchosKokkosComm KokkosKernels
  LIB_OPTIONAL_PACKAGES TpetraTSQR TeuchosNumerics
  LIB_OPTIONAL_TPLS MPI CUDA QD quadmath mpi_advance
)
