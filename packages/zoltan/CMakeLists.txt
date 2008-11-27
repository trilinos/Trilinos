
INCLUDE(TrilinosPackageHelpers)
INCLUDE(Trilinos_Add_Option)

#
# A) Define the package
#

TRILINOS_PACKAGE(Zoltan)

#
# B) Set up package-specific options
#

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_ENABLE_OCT
  HAVE_ZOLTAN_OCT
  "Enable support for OCT in ${PACKAGE_NAME}."
  OFF )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_ENABLE_ParMETIS
  HAVE_PARMETIS
  "Enable support for ParMETIS in ${PACKAGE_NAME}."
  OFF )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_ENABLE_Scotch
  HAVE_SCOTCH
  "Enable support for Scotch in ${PACKAGE_NAME}."
  OFF )

OPTION(${PACKAGE_NAME}_ENABLE_CPPDRIVER
  "Enable support for Scotch in ${PACKAGE_NAME}."
  ${Trilinos_ENABLE_CXX} )

#
# C) Add the libraries, tests, and examples
#

IF (NOT TPL_ENABLE_MPI)
  ADD_SUBDIRECTORY(siMPI)
ENDIF()

ADD_SUBDIRECTORY(src)

TRILINOS_PACKAGE_ADD_TEST_DIRECTORIES(test)

#
# D) Do standard postprocessing
#

TRILINOS_PACKAGE_POSTPROCESS()
