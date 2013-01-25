
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

# ENABLE_SHADOWING_WARNINGS  not possible, because of TBB
TRIBITS_PACKAGE( STK )

#
# B) Set up package-specific options
#

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Trilinos
  HAVE_STK_Trilinos
  "Building as a Trilinos package"
  ON
  )

# set fortran mangling options

IF("${FC_FN_UNDERSCORE}" STREQUAL "NO_UNDER")
  SET(FORTRAN_NO_UNDERSCORE ON)
ELSEIF("${FC_FN_UNDERSCORE}" STREQUAL "UNDER")
  SET(FORTRAN_ONE_UNDERSCORE ON)
ELSEIF("${FC_FN_UNDERSCORE}" STREQUAL "SECOND_UNDER")
  SET(FORTRAN_TWO_UNDERSCORES ON)
ELSE()
  MESSAGE("Could not determine the Fortran mangling defaulting to no underscore.")
  SET(FORTRAN_NO_UNDERSCORE ON)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(stk_util)

ADD_SUBDIRECTORY(stk_mesh)

ADD_SUBDIRECTORY(stk_algsup)

IF (STK_ENABLE_Sacado)
  ADD_SUBDIRECTORY(stk_expreval)
ENDIF()

IF( ${PACKAGE_NAME}_ENABLE_SEACASIoss )
  ADD_SUBDIRECTORY(stk_io)
ENDIF()

IF( STK_ENABLE_FEI AND STK_ENABLE_AztecOO )
  ADD_SUBDIRECTORY(stk_linsys)
ENDIF()

IF( STK_ENABLE_Zoltan AND TPL_ENABLE_MPI )
  ADD_SUBDIRECTORY(stk_rebalance)
ENDIF()

ADD_SUBDIRECTORY(stk_transfer)

IF( STK_ENABLE_ThreadPool )
  ADD_SUBDIRECTORY(stk_search)
  ADD_SUBDIRECTORY(stk_search_util)
ENDIF()

#IF( STK_ENABLE_Intrepid AND STK_ENABLE_SEACASIoss AND STK_ENABLE_BoostLib AND STK_ENABLE_OpenNURBS)
#IF( STK_ENABLE_Intrepid AND STK_ENABLE_SEACASIoss AND STK_ENABLE_OpenNURBS)
IF( STK_ENABLE_Intrepid AND STK_ENABLE_SEACASIoss AND STK_ENABLE_ThreadPool )
  ADD_SUBDIRECTORY(stk_percept)
  ADD_SUBDIRECTORY(stk_adapt)
ENDIF()

IF ( STK_ENABLE_BoostLib AND STK_ENABLE_TESTS )
  ADD_SUBDIRECTORY(stk_usecases)
ENDIF()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
