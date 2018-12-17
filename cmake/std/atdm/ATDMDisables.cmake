#####################################################
###
###  Standrd disables for ATDM Trilinos builds
###
#####################################################

ATDM_SET_CACHE(TPL_ENABLE_GLM OFF CACHE BOOL)
ATDM_SET_CACHE(TPL_ENABLE_Matio OFF CACHE BOOL)
ATDM_SET_CACHE(TPL_ENABLE_SuperLU OFF CACHE BOOL)
ATDM_SET_CACHE(TPL_ENABLE_X11 OFF CACHE BOOL)
ATDM_SET_CACHE(TPL_ENABLE_yaml-cpp OFF CACHE BOOL)

# Packages and sub-packages disables common to both SPARC and EMPIRE
SET(ATDM_SE_PACKAGE_DISABLES
  MiniTensor
  GlobiPack
  OptiPack
  Isorropia
  KokkosExample
  MiniTensor
  TpetraTSQR
  Domi
  Pliris
  Komplex
  Trios
  FEI
  TriKota
  Intrepid
  STKClassic
  STKSearchUtil
  STKUnit_tests
  STKDoc_tests
  STKExp
  Moertel
  ShyLU_DD
  ShyLU
  Stokhos
  MOOCHO
  PyTrilinos
  ForTrilinos
  TrilinosCouplings
  Pike
  )

IF (NOT ATDM_ENABLE_SPARC_SETTINGS)
  # Extra disables for the non-SPARC (EMPIRE) build.
  SET(ATDM_SE_PACKAGE_DISABLES
    ${ATDM_SE_PACKAGE_DISABLES}
    ShyLU_Node
    ROL
    )
  # NOTE: For now, we will disable these packages not being used by EMPIRE for
  # now so that we don't introduce any new failing tests in the existing ATDM
  # Trilinos builds.  But as we get the SPARC ATDM Trilinos configuration
  # building on more machines and we test SPARC and EMPIRE against the fuller
  # SPARC configuration, then ShyLU_Node and ROL will get enabled one machine
  # at a time in an orderly fashion.
ENDIF()


#
# Set of ATDM Trilinos packages for wich we don't want to run the test suite.
#
# This allows us to use Trilinos_ENABLE_ALL_PACKAGES=ON
#
SET(ATDM_SE_PACKAGE_TEST_DISABLES
  TrilinosFrameworkTests
  Gtest
  RTOp
  Epetra
  Zoltan
  Shards
  Triutils
  EpetraExt
  AztecOO
  Galeri
  Amesos
  Pamgen
  Ifpack
  ML
  )

#
# We should have the same set of disables for the SPARC and EMPIRE builds!
# Therefore, let's not add any more disables for now!  This just enables more
# subpackages which should hopefully be harmless.  Also, with the STK disables
# for the SPARC settings, I got build failures in the remaining STK packages.
# We will just need to get SPARC and EMPIRE and Trilinos to work with these
# same sets of enabled packages (even if it is more than they were using
# before).
#

IF (ATDM_ADD_EXTRA_APP_SPECIFIC_DISABLES)  # Undefined and therefore fasle!

IF (ATDM_ENABLE_SPARC_SETTINGS)
  # Add extra disables if using SPARC settings
  SET(ATDM_SE_PACKAGE_DISABLES
    ${ATDM_SE_PACKAGE_DISABLES}
    STKMesh
    STKIO
    STKTopology
    )
ELSE()
  # Add extra disables if not adding SPARC settings
  SET(ATDM_SE_PACKAGE_DISABLES
    ${ATDM_SE_PACKAGE_DISABLES}
    ShyLU_Node
    SEACASExodus_for
    SEACASExoIIv2for32
    SEACASSupes
    SEACASSuplib
    SEACASSVDI
    SEACASPLT
    SEACASBlot
    SEACASConjoin
    SEACASEjoin
    SEACASExo2mat
    SEACASExomatlab
    SEACASExotxt
    SEACASExo_format
    SEACASEx1ex2v2
    SEACASFastq
    SEACASGjoin
    SEACASGen3D
    SEACASGenshell
    SEACASGrepos
    SEACASGrope
    SEACASMapvarlib
    SEACASMapvar
    SEACASMapvar-kd
    SEACASMat2exo
    SEACASNumbers
    SEACASTxtexo
    SEACASEx2ex1v2
    STKUnit_test_utils # This is needed to test STKSearch!
    STKSearch
    STKTransfer
    ROL
    )
ENDIF()

ENDIF (ATDM_ADD_EXTRA_APP_SPECIFIC_DISABLES)

# ToDo: These disables need to be merged into a single list of disables!

# Disable the disabled SE packages
FOREACH(ATDM_SE_PACKAGE ${ATDM_SE_PACKAGE_DISABLES})
  ATDM_SET_CACHE(Trilinos_ENABLE_${ATDM_SE_PACKAGE} OFF CACHE BOOL)
ENDFOREACH()

# Disable the disabled SE package tests and examples
FOREACH(ATDM_SE_PACKAGE ${ATDM_SE_PACKAGE_TEST_DISABLES})
  ATDM_SET_CACHE(${ATDM_SE_PACKAGE}_ENABLE_TESTS OFF CACHE BOOL)
  ATDM_SET_CACHE(${ATDM_SE_PACKAGE}_ENABLE_EXAMPLES OFF CACHE BOOL)
ENDFOREACH()
