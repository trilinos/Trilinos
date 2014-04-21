SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  #SubPackageName       Directory         Class    Req/Opt
  #
  # New Kokkos subpackages:
  TPL                   TPL               PS       OPTIONAL
  Core                  core              SS       OPTIONAL
  Compat                compat            SS       OPTIONAL
  Classic               classic           PS       OPTIONAL
  Containers            containers        SS       OPTIONAL
  LinAlg                linalg            SS       OPTIONAL
  Example               example           SS       OPTIONAL
  MpiComm               mpicomm           SS       OPTIONAL
  Task                  task              EX       OPTIONAL
  )

SET(LIB_REQUIRED_DEP_PACKAGES )
SET(LIB_OPTIONAL_DEP_PACKAGES )
SET(TEST_REQUIRED_DEP_PACKAGES )
SET(TEST_OPTIONAL_DEP_PACKAGES )
SET(LIB_REQUIRED_DEP_TPLS )
SET(LIB_OPTIONAL_DEP_TPLS ) 
SET(TEST_REQUIRED_DEP_TPLS )
SET(TEST_OPTIONAL_DEP_TPLS )
