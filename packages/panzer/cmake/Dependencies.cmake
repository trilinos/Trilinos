SET(LIB_REQUIRED_DEP_PACKAGES Teuchos Kokkos Sacado Phalanx Intrepid ThyraCore ThyraTpetraAdapters ThyraEpetraAdapters Tpetra Epetra EpetraExt Zoltan Stratimikos Piro NOX Rythmos KokkosCore)
SET(LIB_OPTIONAL_DEP_PACKAGES STKClassic SEACASIoss SEACASExodus Teko MueLu Ifpack2 FEI)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES Stratimikos Piro NOX Rythmos)
SET(LIB_REQUIRED_DEP_TPLS MPI BoostLib)
SET(LIB_OPTIONAL_DEP_TPLS PAPI)

# The Boostlib dependency below should not be there for tests as it is
# already in the lib requirement above, but there is an issue in
# tribits for when treating TPLs like system libraries for warning as
# error builds.  The addition of BoostLib below is a hack to suppress
# boost warnings from warning as error builds in tests.  When TriBITS
# issue 63 is addressed, this issue will go away.

SET(TEST_REQUIRED_DEP_TPLS BoostLib)
SET(TEST_OPTIONAL_DEP_TPLS)

TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(MueLu)
