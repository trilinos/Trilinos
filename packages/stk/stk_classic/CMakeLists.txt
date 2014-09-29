

#
# A) Define the package
#

# ENABLE_SHADOWING_WARNINGS  not possible, because of TBB
TRIBITS_SUBPACKAGE(Classic)

#
# B) Set up package-specific options
#

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Trilinos
  HAVE_STK_Trilinos
  "Building as a Trilinos package"
  ON
  )

## set fortran mangling options
#
#IF("${FC_FN_UNDERSCORE}" STREQUAL "NO_UNDER")
#  SET(FORTRAN_NO_UNDERSCORE ON)
#ELSEIF("${FC_FN_UNDERSCORE}" STREQUAL "UNDER")
#  SET(FORTRAN_ONE_UNDERSCORE ON)
#ELSEIF("${FC_FN_UNDERSCORE}" STREQUAL "SECOND_UNDER")
#  SET(FORTRAN_TWO_UNDERSCORES ON)
#ELSE()
#  MESSAGE("Could not determine the Fortran mangling defaulting to no underscore.")
#  SET(FORTRAN_NO_UNDERSCORE ON)
#ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(stk_util)

ADD_SUBDIRECTORY(stk_mesh)

ADD_SUBDIRECTORY(stk_algsup)

IF (${PACKAGE_NAME}_ENABLE_Sacado)
  ADD_SUBDIRECTORY(stk_expreval)
ENDIF()

IF( ${PACKAGE_NAME}_ENABLE_SEACASIoss )
  ADD_SUBDIRECTORY(stk_io)
ENDIF()

IF( ${PACKAGE_NAME}_ENABLE_FEI AND ${PACKAGE_NAME}_ENABLE_AztecOO )
  ADD_SUBDIRECTORY(stk_linsys)
ENDIF()

IF( ${PACKAGE_NAME}_ENABLE_Zoltan AND TPL_ENABLE_MPI )
  ADD_SUBDIRECTORY(stk_rebalance)
ENDIF()

ADD_SUBDIRECTORY(stk_transfer)

IF( ${PACKAGE_NAME}_ENABLE_ThreadPool )
  ADD_SUBDIRECTORY(stk_search)
  ADD_SUBDIRECTORY(stk_search_util)
ENDIF()

IF( ${PACKAGE_NAME}_ENABLE_Intrepid AND ${PACKAGE_NAME}_ENABLE_SEACASIoss AND ${PACKAGE_NAME}_ENABLE_ThreadPool )
  ADD_SUBDIRECTORY(stk_percept)
  ADD_SUBDIRECTORY(stk_adapt)
ENDIF()

IF ( ${PACKAGE_NAME}_ENABLE_BoostLib AND ${PACKAGE_NAME}_ENABLE_TESTS )
  ADD_SUBDIRECTORY(stk_usecases)
ENDIF()

#
# D) Do standard postprocessing
#

TRIBITS_SUBPACKAGE_POSTPROCESS()
