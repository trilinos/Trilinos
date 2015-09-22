SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  #SubPackageName       Directory         Class    Req/Opt
  #
  # New Intrepid2
  Core                  core               PS      REQUIRED
  Intrepid2             intrepid2          EX      OPTIONAL
  )
SET(LIB_REQUIRED_DEP_PACKAGES Teuchos Shards Sacado)
SET(LIB_OPTIONAL_DEP_PACKAGES KokkosCore KokkosAlgorithms)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES Epetra EpetraExt Amesos Pamgen)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS Boost)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
