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

# Packages and sub-packages disables.
SET(ATDM_SE_PACKAGE_DISABLES
  ThreadPool
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
# NOTE: The above list came from the "Final set of non-enabled SE packages"
# from the output for using the
# em-plasma/BuildScripts/ALL/configure-trilinos.sh script (see TRIL-171).
# This list can be easily maintained going forward.

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

# ToDo: These disables need to be merged into a single list of disables!

# Disable the disabled SE packages
FOREACH(ATDM_SE_PACKAGE ${ATDM_SE_PACKAGE_DISABLES})
  ATDM_SET_CACHE(Trilinos_ENABLE_${ATDM_SE_PACKAGE} OFF CACHE BOOL)
ENDFOREACH()
