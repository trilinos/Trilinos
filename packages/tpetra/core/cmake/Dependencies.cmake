TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_PACKAGES Teuchos TpetraClassic KokkosCore KokkosContainers KokkosAlgorithms TeuchosKokkosCompat TeuchosKokkosComm KokkosKernels
  LIB_OPTIONAL_PACKAGES Epetra TpetraTSQR TeuchosNumerics
  LIB_OPTIONAL_TPLS MPI CUDA QD quadmath
)

IF(TPL_ENABLE_CUDA)
  TRIBITS_TPL_TENTATIVELY_ENABLE(CUSPARSE)
ENDIF()
