TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_PACKAGES Teuchos TpetraClassic
  LIB_OPTIONAL_PACKAGES Epetra KokkosCore KokkosCompat KokkosLinAlg KokkosContainers KokkosAlgorithms KokkosMpiComm TpetraTSQR
  LIB_OPTIONAL_TPLS MPI CUDA Thrust QD Cusp
)
