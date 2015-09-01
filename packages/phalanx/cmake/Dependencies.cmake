SET(LIB_REQUIRED_DEP_PACKAGES Teuchos KokkosCore Sacado)
SET(LIB_OPTIONAL_DEP_PACKAGES Intrepid)
SET(TEST_REQUIRED_DEP_PACKAGES Shards)
#SET(TEST_OPTIONAL_DEP_PACKAGES Belos Epetra Ifpack ML Amesos)
SET(TEST_OPTIONAL_DEP_PACKAGES )
SET(LIB_REQUIRED_DEP_TPLS Boost)
SET(LIB_OPTIONAL_DEP_TPLS )

# The Boost dependency below should not be there for tests as it is
# already in the lib requirement above, but there is an issue in
# tribits for when treating TPLs like system libraries for
# warning-as-error builds.  The addition of Boost below is a temporary
# hack to suppress boost warnings from warning-as-error builds in
# tests.  When TriBITS issue 63 is addressed, this issue will go away.

SET(TEST_REQUIRED_DEP_TPLS Boost)
SET(TEST_OPTIONAL_DEP_TPLS TVMET)
