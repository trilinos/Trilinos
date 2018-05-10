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
  MiniTensor
  GlobiPack
  OptiPack
  Isorropia
  KokkosExample
  MiniTensor
  TpetraTSQR
  Isorropia
  ShyLU_Node
  SEACASExodus_for
  SEACASExoIIv2for32
  SEACASSupes
  SEACASSuplib
  SEACASSVDI
  SEACASPLT
  SEACASAlgebra
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
  Trios
  FEI
  TriKota
  Intrepid
  STKClassic
  STKUnit_test_utils
  STKSearch
  STKSearchUtil
  STKTransfer
  STKUnit_tests
  STKDoc_tests
  STKExp
  Moertel
  ShyLU_DD
  ShyLU
  Stokhos
  MOOCHO
  ROL
  PyTrilinos
  ForTrilinos
  )
# NOTE: The above list came from the "Final set of non-enabled SE packages"
# from the output for using the
# em-plasma/BuildScripts/ALL/configure-trilinos.sh script (see TRIL-171).
# This list can be easily maintained going forward.

# Disable the disabled SE packages
FOREACH(ATDM_SE_PACKAGE ${ATDM_SE_PACKAGE_DISABLES})
  ATDM_SET_CACHE(Trilinos_ENABLE_${ATDM_SE_PACKAGE} OFF CACHE BOOL)
ENDFOREACH()
