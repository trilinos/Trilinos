
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Stokhos)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  STOKHOS_DEBUG
  "Enable debug code in sacado"
  ${Trilinos_ENABLE_DEBUG} )

#IF(Stokhos_ENABLE_ML)
#  SET(HAVE_STOKHOS_ML ON)
#ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
