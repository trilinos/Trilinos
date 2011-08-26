#
# We list MueLu as a dependency, but we only use Xpetra.  So
#   when Xpetra is split out as its own package, we will
#   depend on Xpetra only.
#
SET(LIB_REQUIRED_DEP_PACKAGES Tpetra Kokkos Epetra Teuchos MueLu)
SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_REQUIRED_DEP_PACKAGES Tpetra Kokkos Epetra Teuchos)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS METIS ParMETIS Scotch CCOLAMD)
SET(TEST_REQUIRED_DEP_TPLS) 
SET(TEST_OPTIONAL_DEP_TPLS METIS ParMETIS Scotch CCOLAMD)

