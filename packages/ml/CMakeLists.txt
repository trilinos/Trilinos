
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)
INCLUDE(DualScopeSet)

#
# A) Define the package
#

PACKAGE(ML)

#
# B) Set up package-specific options
#

IF (${PACKAGE_NAME}_ENABLE_Zoltan AND ${PACKAGE_NAME}_ENABLE_ZoltanTpl)
  MESSAGE( FATAL_ERROR
        "\n***"
        "\n*** ERROR:  You can not enable both the package Zoltan and the TPL ZoltanTpl at the same time!"
        "\n***\n"
        )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Zoltan AND NOT ${PACKAGE_NAME}_ENABLE_MPI)
  MESSAGE(STATUS "  WARNING:  Setting ${PACKAGE_NAME}_ENABLE_Zoltan=OFF since ${PACKAGE_NAME}_ENABLE_MPI=OFF!")
  DUAL_SCOPE_SET(${PACKAGE_NAME}_ENABLE_Zoltan "")
  DUAL_SCOPE_SET(HAVE_ML_ZOLTAN OFF)
ENDIF()

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Aztec
  HAVE_ML_AZTEC
  "Enable the Aztec functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Aztec2_1
  HAVE_ML_AZTEC2_1
  "Enable the Aztec 2.1 functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SUPERLU
  HAVE_ML_SUPERLU
  "Enable the SUPERLU functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SUPERLU2_0
  HAVE_ML_SUPERLU2_0
  "Enable the SUPERLU 2.0 functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SUPERLUDIST
  HAVE_ML_SUPERLUDIST
  "Enable the SUPERLUDIST functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Enrich
  HAVE_ML_ENRICH
  "Enable the Enrich functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Memory_Check
  HAVE_ML_MEMORY_CHECK
  "Enable the memory checking functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_New_T_PE
  HAVE_ML_NEW_T_PE
  "Enable the new T PE functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Complex_Maxwell
  HAVE_ML_COMPLEX_MAXWELL
  "Enable the complex Maxwell functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Timing
  HAVE_ML_TIMING
  "Enable the timing functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Flops
  HAVE_ML_FLOPS
  "Enable the flops functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Flops
  HAVE_ML_FLOPS
  "Enable the flops functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_MLapi
  HAVE_ML_MLAPI
  "Enable the ML api."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Cfunc
  HAVE_ML_CFUNC
  "Enable the Cfunc functionality."
  OFF )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(examples)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
