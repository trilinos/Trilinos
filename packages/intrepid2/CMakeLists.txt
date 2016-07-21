

#
# A) Define the package
#

TRIBITS_PACKAGE(Intrepid2)

#
# B) Set up package-specific options
#

IF(Amesos_ENABLE_UMFPACK)
  SET(${PACKAGE_NAME}_ENABLE_AMESOS_UMFPACK_DEFAULT ON)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_AMESOS_UMFPACK_DEFAULT OFF)
ENDIF()

IF(Trilinos_ENABLE_ROL)
  SET(${PACKAGE_NAME}_ENABLE_ROL ON)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_ROL OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_AMESOS_UMFPACK
  HAVE_INTREPID2_AMESOS_UMFPACK
  "Enable Amesos with UmfPack functionality."
  ${Intrepid2_ENABLE_AMESOS_UMFPACK_DEFAULT} )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  HAVE_INTREPID2_DEBUG
  "Enable a host of runtime debug checking."
  ${${PROJECT_NAME}_ENABLE_DEBUG}
  )

IF(${PACKAGE_NAME}_ENABLE_DEBUG_INF_CHECK_OVERRIDE)
  ## override the inf flag when it is needed e.g., refactoring dyn rank view
  TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG_INF_CHECK
    HAVE_INTREPID2_DEBUG_INF_CHECK
    "Enable a host of runtime debug checking for infs and nans."
    ${${PROJECT_NAME}_ENABLE_DEBUG})
ELSE()
  ## unless override is enabled, inf check is turned off by default
  TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG_INF_CHECK
    HAVE_INTREPID2_DEBUG_INF_CHECK
    "Enable a host of runtime debug checking for infs and nans."
    OFF)
ENDIF()   

#
# C) Add the libraries, tests, and examples
#

IF (Intrepid2_ENABLE_KokkosDynRankView)
ADD_SUBDIRECTORY(refactor)
ELSE()
ADD_SUBDIRECTORY(core)
ENDIF()

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
