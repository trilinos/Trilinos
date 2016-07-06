SET(LIB_REQUIRED_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_PACKAGES
  EpetraExt Isorropia Amesos AztecOO Belos Ifpack ML NOX Zoltan STKClassic STK Stokhos)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES Amesos Isorropia Epetra EpetraExt Ifpack Intrepid Pamgen AztecOO ML Zoltan STKClassic Teko TpetraKernels Tpetra MueLu KokkosCore TeuchosKokkosCompat KokkosContainers TeuchosKokkosComm Stokhos)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS TASMANIAN)

# STKClassic support in TrilinosCouplings is broken (see Trilinos #19)
SET(TrilinosCouplings_ENABLE_STKClassic OFF CACHE BOOL
  "Set by default in TrilinosCopulings Dependencies.cmake")
