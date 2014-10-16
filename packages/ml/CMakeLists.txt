
INCLUDE(DualScopeSet)

#
# A) Define the package
#

TRIBITS_PACKAGE(ML  DISABLE_CIRCULAR_REF_DETECTION_FAILURE)

#
# B) Set up package-specific options
#

IF (${PACKAGE_NAME}_ENABLE_Zoltan AND NOT ${PACKAGE_NAME}_ENABLE_MPI)
  MESSAGE(STATUS "  WARNING:  Setting ${PACKAGE_NAME}_ENABLE_Zoltan=OFF since ${PACKAGE_NAME}_ENABLE_MPI=OFF!")
  DUAL_SCOPE_SET(${PACKAGE_NAME}_ENABLE_Zoltan OFF)
  DUAL_SCOPE_SET(HAVE_ML_ZOLTAN OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Aztec
  HAVE_ML_AZTEC
  "Enable the Aztec functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SuperLU1_0
  HAVE_ML_SUPERLU
  "Enable the SUPERLU functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SuperLU2_0
  HAVE_ML_SUPERLU2_0
  "Enable the SUPERLU 2.0 functionality."
  OFF )

# macro ADD_OPTION_AND_DEFINE should not be used here:
# because SuperLU is a TPL the CMake build system defines ${PACKAGE_NAME}_ENABLE_SuperLU
IF(${PACKAGE_NAME}_ENABLE_SuperLU) #TODO: can it be safely replaced by: SET(HAVE_ML_SUPERLU4_0 "${${PACKAGE_NAME}_ENABLE_SuperLU}") ?
   SET(HAVE_ML_SUPERLU4_0 ON)
ELSE()
   SET(HAVE_ML_SUPERLU4_0 OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SUPERLUDIST
  HAVE_ML_SUPERLUDIST
  "Enable the SUPERLUDIST functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Enrich
  HAVE_ML_ENRICH
  "Enable the Enrich functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Memory_Check
  HAVE_ML_MEMORY_CHECK
  "Enable the memory checking functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_New_T_PE
  HAVE_ML_NEW_T_PE
  "Enable the new T PE functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Complex_Maxwell
  HAVE_ML_COMPLEX_MAXWELL
  "Enable the complex Maxwell functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Timing
  HAVE_ML_TIMING
  "Enable the timing functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Flops
  HAVE_ML_FLOPS
  "Enable the flops functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_MLapi
  HAVE_ML_MLAPI
  "Enable the ML api."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Cfunc
  HAVE_ML_CFUNC
  "Enable the Cfunc functionality."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TekoSmoothers
  HAVE_ML_TekoSmoothers
  "Enable Teko Smoothing Functionality."
  OFF )

# Note: macro HAVE_ML_MATLAB is automatically defined if TPL MATLAB is enabled.

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(matlab)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(examples)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  util
  etc
  test/tmp
  doc/UsersGuide
  doc/DevelopersGuide
  doc/MLAPI
  python
  doc/DoxyfileWeb
  doc/build_docs
  doc/ml-logo.eps
  doc/ml-logo.jpg
  doc/sc2000.ps.gz
  examples/Makefile-common.include
  examples/Maxwell/ml_periodic_max.c
  examples/Other/ml_read_complexmaxwell.c
  examples/Other/ml_read_maxwell.c
  examples/Other/ml_star2d.c
  examples/Other/new_readex.c
  examples/Other/oldml_readex.c
  examples/Other/seg_readex.c
  examples/README.AddingExamples
  examples/RefMaxwell
  examples/RefMaxwell/rpc.cpp
  src/Coarsen/README
  src/Main/ml_v_cycle.c
  src/Smoother/README
  src/Utils/jmpilib.c
  src/Utils/jostle.h
  src/Utils/ml_vampir.c
  src/Utils/tumi.c
  test/README.runtests
  test/Zoltan/cxx_main_simple.cpp
  test/scripts
  test/scripts/run-tests.sh
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
