
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Zoltan2)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_METIS
  HAVE_METIS
  "Metis is present"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_ParMETIS
  HAVE_PARMETIS
  "ParMetis is present"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Scotch
  HAVE_SCOTCH
  "Scotch is present"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_PTScotch
  HAVE_PTSCOTCH
  "PTScotch is present"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_NEMESIS
  HAVE_NEMESIS
  "Nemesis/Exodus is present"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_HUND
  HAVE_HUND
  "HUND is present"
  OFF )

# 
# Include Fortran support and Fortran tests
#
# IF(Trilinos_ENABLE_Fortran)
# ENDIF()
#

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
