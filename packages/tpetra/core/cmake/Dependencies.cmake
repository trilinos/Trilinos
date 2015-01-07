TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_PACKAGES Teuchos TpetraClassic
  LIB_OPTIONAL_PACKAGES Epetra KokkosCore KokkosContainers KokkosAlgorithms TeuchosKokkosCompat TeuchosKokkosComm TpetraKernels TpetraTSQR
  LIB_OPTIONAL_TPLS MPI CUDA Thrust QD Cusp
)
