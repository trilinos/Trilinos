
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Intrepid)

#
# B) Set up package-specific options
#

IF(Amesos_ENABLE_UMFPACK)
  SET(${PACKAGE_NAME}_ENABLE_AMESOS_UMFPACK_DEFAULT ON)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_AMESOS_UMFPACK_DEFAULT OFF)
ENDIF()

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_AMESOS_UMFPACK
  HAVE_INTREPID_AMESOS_UMFPACK
  "Enable Amesos with UmfPack functionality."
  ${Intrepid_ENABLE_AMESOS_UMFPACK_DEFAULT} )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  HAVE_INTREPID_DEBUG
  "Enable a host of runtime debug checking."
  ${Trilinos_ENABLE_DEBUG}
  )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG_INF_CHECK
  HAVE_INTREPID_DEBUG_INF_CHECK
  "Enable a host of runtime debug checking for infs and nans."
  ${Trilinos_ENABLE_DEBUG}
  )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
