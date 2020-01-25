#####################################################
###
###  Standrd disables for ATDM Trilinos builds
###
#####################################################


#
# A) TPL disables
#

ATDM_SET_CACHE(TPL_ENABLE_GLM OFF CACHE BOOL)
ATDM_SET_CACHE(TPL_ENABLE_Matio OFF CACHE BOOL)
ATDM_SET_CACHE(TPL_ENABLE_SuperLU OFF CACHE BOOL)
ATDM_SET_CACHE(TPL_ENABLE_X11 OFF CACHE BOOL)
ATDM_SET_CACHE(TPL_ENABLE_yaml-cpp OFF CACHE BOOL)


#
# B) SE Package disables
#
# These are SE packages that we don't even want enabled.  (This is the
# "black-listing" approach)
#

# Packages and sub-packages disables common to both SPARC and EMPIRE
SET(ATDM_SE_PACKAGE_DISABLES
  TrilinosFrameworkTests
  MiniTensor
  Isorropia
  KokkosExample
  Domi
  Pliris
  Komplex
  FEI
  TriKota
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

IF (ATDM_USE_CUDA AND ATDM_COMPLEX)

  # Disable extra packages not being used by GEMMA in cuda+complex builds (see
  # ATDV-263, ATDV-265)
  SET(ATDM_SE_PACKAGE_DISABLES
    ${ATDM_SE_PACKAGE_DISABLES}
    ShyLU_Node
    Amesos2
    SEACAS
    Anasazi
    Ifpack2
    Stratimikos
    Teko
    Intrepid
    Intrepid2
    STK
    Percept
    NOX
    Moertel
    MueLu
    Rythmos
    Tempus
    ROL
    Piro
    Panzer
    )
  # Note, above:
  # * We are allowing the enable of Zoltan2 for now even though GEMMA is not
  #   currently using it because Zoltan2 developers want to maintain it (see
  #   discussion in ATDV-263).
  # * The package MueLu does not currently build for cuda+complex builds so
  #   therefore must be disabled (see #4599)

ENDIF()


#
# C) SE Package test disables
#
# This is the is the ATDM Trilinos SE packages that we allow to be enabled but
# for wich we don't want to run the test suite.
#
# This allows us to use Trilinos_ENABLE_ALL_PACKAGES=ON to enable the ATDM
# Trilinos builds and get the tests for only the packages we want (i.e the
# "black-listing" appraoch.
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
  Anasazi
  Intrepid
  Piro
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


#
# D) Individual Executable and Test Disables
#
# There are some package tests that have to be disabled for a broad set of
# builds for example, if all openmp builds are failing a certain test then it
# makes more sense to disbale it once in this file instead of in every openmp
# buid's tweaks file
#

# Issue #3638
ATDM_SET_ENABLE(Teko_ModALPreconditioner_MPI_1_DISABLE ON)

# Disable Zoltan2_XpetraEpertraMatrix exec that does not build with no global
# int instatiation (see #5411)
ATDM_SET_ENABLE(Zoltan2_XpetraEpetraMatrix_EXE_DISABLE ON)
ATDM_SET_ENABLE(Zoltan2_XpetraEpetraMatrix_MPI_4_DISABLE ON)

# Diable experimental Belos HybridGMES tests in all ATDM Trilinos builds as
# none of the ATDM customer codes are using this solver (see #4159)
ATDM_SET_ENABLE(Belos_Tpetra_HybridGMRES_hb_test_1_MPI_4_DISABLE ON)
ATDM_SET_ENABLE(Belos_Tpetra_HybridGMRES_hb_test_0_MPI_4_DISABLE ON)

# Disable test for Belos recycling CG solver that is not used by ATDM APPs
# (see #2919)
ATDM_SET_ENABLE(Belos_rcg_hb_MPI_4_DISABLE ON)

# Disable Piro_ThyraSolver exec that does not build with no global int
# instantiation (see #5412)
ATDM_SET_ENABLE(Piro_ThyraSolver_EXE_DISABLE ON)
ATDM_SET_ENABLE(Piro_ThyraSolver_MPI_4_DISABLE ON)

# Disable Piro_AnalysisDriverTpetra exec that will not build with no global
# int instantiation (see #5446)
ATDM_SET_ENABLE(Piro_AnalysisDriverTpetra_EXE_DISABLE ON)
ATDM_SET_ENABLE(Piro_AnalysisDriverTpetra_MPI_4_DISABLE ON)

# Disable ROL test exec that will not buld with no global int instantiation
# (see #5447)
ATDM_SET_ENABLE(ROL_adapters_tpetra_test_vector_SimulatedVectorTpetraBatchManagerInterface_EXE_DISABLE ON)
ATDM_SET_ENABLE(ROL_adapters_tpetra_test_vector_SimulatedVectorTpetraBatchManagerInterface_MPI_4_DISABLE ON)

IF ("${ATDM_CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")

  ATDM_SET_ENABLE(PanzerAdaptersSTK_CurlLaplacianExample-ConvTest-Quad-Order-4_DISABLE ON)
  ATDM_SET_ENABLE(PanzerAdaptersSTK_MixedPoissonExample-ConvTest-Hex-Order-3_DISABLE ON)

ENDIF()

IF (ATDM_NODE_TYPE STREQUAL "CUDA")

  # Disable ROL tests that don't work with CUDA builds (see #3543, #6124)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_0ld_adv-diff-react_example_01_MPI_4_DISABLE ON)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_0ld_adv-diff-react_example_02_MPI_4_DISABLE ON)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_0ld_poisson_example_01_MPI_4_DISABLE ON)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_0ld_stefan-boltzmann_example_03_MPI_4_DISABLE ON)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_navier-stokes_example_01_MPI_4_DISABLE ON)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_navier-stokes_example_02_MPI_4_DISABLE ON)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_nonlinear-elliptic_example_01_MPI_4_DISABLE ON)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_nonlinear-elliptic_example_02_MPI_4_DISABLE ON)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_obstacle_example_01_MPI_4_DISABLE ON)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_stefan-boltzmann_example_01_MPI_4_DISABLE ON)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_stefan-boltzmann_example_03_MPI_4_DISABLE ON)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_topo-opt_poisson_example_01_MPI_4_DISABLE ON)
  ATDM_SET_ENABLE(ROL_test_elementwise_TpetraMultiVector_MPI_4_DISABLE ON)

  # Disable Zoltan tests (see #3749)
  ATDM_SET_ENABLE(TrilinosCouplings_Example_Maxwell_MueLu_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(TrilinosCouplings_Example_Maxwell_MueLu_MPI_4_DISABLE ON)

  # Disable Zoltan tests (see #4042)
  ATDM_SET_ENABLE(Zoltan_ch_ewgt_zoltan_parallel_DISABLE ON)
  ATDM_SET_ENABLE(Zoltan_ch_grid20x19_zoltan_parallel_DISABLE ON)
  ATDM_SET_ENABLE(Zoltan_ch_nograph_zoltan_parallel_DISABLE ON)
  ATDM_SET_ENABLE(Zoltan_ch_simple_zoltan_parallel_DISABLE ON)

ENDIF()
