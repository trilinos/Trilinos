
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Zoltan)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_OCT
  HAVE_ZOLTAN_OCT
  "Enable support for OCT in ${PACKAGE_NAME}."
  ON )

OPTION(${PACKAGE_NAME}_ENABLE_CPPDRIVER
  "Enable C++ driver for ${PACKAGE_NAME}."
  ${Trilinos_ENABLE_CXX} )

#
# C) Add the libraries, tests, and examples
#

IF (NOT TPL_ENABLE_MPI)
  ADD_SUBDIRECTORY(siMPI)
ENDIF()

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(src/driver)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
