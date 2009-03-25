
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Amesos)

#
# B) Set up package-specific options
#

# if using SuperLUDist, must also link in ParMETIS for some reason
IF(${PACKAGE_NAME}_ENABLE_SuperLUDist AND NOT ${PACKAGE_NAME}_ENABLE_ParMETIS)
  MESSAGE(FATAL_ERROR "The Amesos support for the SuperLUIDist TPL requires the ParMETIS TPL.  Either disable Amesos SuperLUDist support or enable the ParMETIS TPL.")
ENDIF()

## Most of the options below are bogus.  The handling of TPLs is done
## through the dependencies file.  This needs to be cleaned up. I'm
## leaving them for now so that we have a list of all options.  They
## will be removed as correct TPL support is added.

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SUPERLU
  HAVE_AMESOS_SUPERLU
  "Enable SuperLU functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DSCPACK
  HAVE_AMESOS_DSCPACK
  "Enable DSCPACK functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_UMFPACK
  HAVE_AMESOS_UMFPACK
  "Enable UMFPACK functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_MC64
  HAVE_AMESOS_MC64
  "Enable MC64 functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_KLU
  HAVE_AMESOS_KLU
  "Enable KLU (Tim Davis's Files)."
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Mumps
  HAVE_AMESOS_MUMPS
  "Enable MUMPS functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SCALAPACK
  HAVE_AMESOS_SCALAPACK
  "Enable SCALAPACK functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Pardiso
  HAVE_AMESOS_PARDISO
  "Enable PARDISO functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TAUCS
  HAVE_AMESOS_TAUCS
  "Enable TAUCS functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_PARAKLETE
  HAVE_AMESOS_PARAKLETE
  "Enable PARAKLETE functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_LAPACK
  HAVE_AMESOS_LAPACK
  "Enable LAPACK functionality."
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Galeri
  HAVE_AMESOS_Galeri
  "Enable Galeri functionality for tests."
  ON )

# ...

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
