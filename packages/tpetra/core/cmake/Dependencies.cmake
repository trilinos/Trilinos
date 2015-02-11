TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_PACKAGES Teuchos TpetraClassic KokkosCore KokkosContainers KokkosAlgorithms TeuchosKokkosCompat TeuchosKokkosComm TpetraKernels
  LIB_OPTIONAL_PACKAGES Epetra TpetraTSQR
  LIB_OPTIONAL_TPLS MPI CUDA QD
)
