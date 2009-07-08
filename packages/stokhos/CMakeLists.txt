
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

SET(Stokhos_ENABLE_Quadrature_Default FALSE)
IF(TPL_ENABLE_Dakota OR TPL_ENABLE_ForUQTK)
  SET(Stokhos_ENABLE_Quadrature_Default TRUE)
ENDIF()
SET(Stokhos_ENABLE_Quadrature ${Stokhos_ENABLE_Quadrature_Default} CACHE BOOL
    "Whether to enable quadrature routines in Stokhos")

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
