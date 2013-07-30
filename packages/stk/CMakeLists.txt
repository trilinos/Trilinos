
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

IF(${FC_FN_UNDERSCORE} STREQUAL "NO_UNDER")
  SET(FORTRAN_NO_UNDERSCORE ON)
ELSEIF(${FC_FN_UNDERSCORE} STREQUAL "UNDER")
  SET(FORTRAN_ONE_UNDERSCORE ON)
ELSEIF(${FC_FN_UNDERSCORE} STREQUAL "SECOND_UNDER")
  SET(FORTRAN_TWO_UNDERSCORES ON)
ELSE()
  MESSAGE(FATAL_ERROR "Could not determine the Fortran mangling.")
ENDIF()

#
# C) Add the libraries, tests, and examples
#

#################################################################################
# NOTE: to build stk_2.0 (new Sierra forked version of stk, forked in
# December 2012), create a symbolic link in this directory to your
# Sierra version of stk.  For example, 'ln -s /scratch/user/code/stk stk_2.0'
# Next, uncomment the next line and comment the one following it.
# Then point your Trilinos configuration to this installation of
# Trilinos, as usual, and configure and build.
# Note the local .gitignore ignores the symlink.
#################################################################################
SET(STK_IS_V2 "ON")
#SET(STK_IS_V2 "OFF")

IF(${STK_IS_V2} STREQUAL "ON")
  MESSAGE("STK_IS_V2 is set")
  SET(STK_ALT_DIR "stk_2.0/")
  ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_topology)
  # FIXME - need to remove this and get the code running without it, to export to Trilinos/open-source 
  add_definitions(-DSIERRA_MIGRATION)

ELSE()
  SET(STK_ALT_DIR "")
ENDIF()

ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_util)

ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_mesh)

ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_algsup)

IF (STK_ENABLE_Sacado)
  ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_expreval)
ENDIF()

IF( ${PACKAGE_NAME}_ENABLE_SEACASIoss )
  ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_io)
ENDIF()

IF( STK_ENABLE_FEI AND STK_ENABLE_AztecOO )
  ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_linsys)
ENDIF()

IF( STK_ENABLE_Zoltan AND TPL_ENABLE_MPI )
  ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_rebalance)
ENDIF()

ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_transfer)

IF( STK_ENABLE_ThreadPool )
  ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_search)
  ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_search_util)
ENDIF()

#IF( STK_ENABLE_Intrepid AND STK_ENABLE_SEACASIoss AND STK_ENABLE_BoostLib AND STK_ENABLE_OpenNURBS)
#IF( STK_ENABLE_Intrepid AND STK_ENABLE_SEACASIoss AND STK_ENABLE_OpenNURBS)
IF( STK_ENABLE_Intrepid AND STK_ENABLE_SEACASIoss AND STK_ENABLE_ThreadPool )
  ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_percept)
  ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_adapt)
ENDIF()

IF ( STK_ENABLE_BoostLib AND STK_ENABLE_TESTS )
  ADD_SUBDIRECTORY(${STK_ALT_DIR}stk_usecases)
ENDIF()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
