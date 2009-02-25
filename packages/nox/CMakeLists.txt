INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)
INCLUDE(CombinedOption)

#
# A) Define the package
#

PACKAGE(NOX)

#
# B) Set up package-specific options
#

# LOCA support
OPTION(${PACKAGE_NAME}_ENABLE_LOCA
  "Enable parameter LOCA support"
  ON )

# MF support
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_MF
  HAVE_LOCA_MF
  "Enable Multifario (MF) support in LOCA"
  OFF )

# Ignoring these flags - should be fixed in more recent compilers
# AC_CXX_STL_MAP_KEY_CONST
# AC_CXX_INVALID_TEMPLATE_QUALIFIER

###################################################
# isnan and isinf support
###################################################
# Now moved into Trilinos top level environment probes

###################################################
# Check for DOT program for Doxygen
###################################################
SET(HAVE_DOXYGEN_DOT NO)
IF(DOXYGEN_DOT_EXECUTABLE)
  SET(HAVE_DOXYGEN_DOT YES)
ENDIF()
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in
               ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile @ONLY)

###################################################
# Debug build
###################################################
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  HAVE_NOX_DEBUG
  "Build with debug code enabled."
  ${Trilinos_ENABLE_DEBUG} )

###################################################
# Optional Trilinos Library Support
###################################################
IF(NOX_ENABLE_ML)
  SET(HAVE_NOX_ML_EPETRA ON)
ENDIF()

IF(NOX_ENABLE_EpetraExt)
  SET(HAVE_NOX_EPETRAEXT ON)
ENDIF()

IF(NOX_ENABLE_Amesos)
  SET(HAVE_NOX_AMESOS ON)
ENDIF()

IF(NOX_ENABLE_Sacado)
  SET(HAVE_NOX_SACADO ON)
ENDIF()

IF(NOX_ENABLE_Anasazi)
  SET(HAVE_LOCA_ANASAZI ON)
ENDIF()

###################################################
# Optional TPL Support
###################################################

###################################################
# Concrete implementations of NOX::Abstract Classes
###################################################

# LAPACK
OPTION(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK
  "Enables support for LAPACK implementation of Abstract Group and Vector classes.  Defaults to OFF"
  OFF
  )

# Epetra
COMBINED_OPTION( NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA
  DEP_OPTIONS_NAMES NOX_ENABLE_AztecOO NOX_ENABLE_Epetra NOX_ENABLE_Ifpack
  DOCSTR "Enables support for Epetra implementation of Abstract Group and"
    " Vector classes.  Defaults to ON if Epetra, AztecOO, and Ifpack are"
    " enabled, otherwise defaults to OFF."
  )

# Thyra
COMBINED_OPTION( NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA
  DEP_OPTIONS_NAMES NOX_ENABLE_Thyra
  DOCSTR "Enables support for Thyra implementation of Abstract Group"
    " and Vector classes.  Defaults to ON if Thyra is enabled,"
    " otherwise defaults to OFF"
  )

###################################################
###################################################

#
# C) Add the libraries, tests, and examples
#

# Limit the number of MPI procs to 2 for now (see bug 4335)
#SET( MPI_EXEC_MAX_NUMPROCS 2 )

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(src-lapack)
ADD_SUBDIRECTORY(src-epetra)
#ADD_SUBDIRECTORY(src-belos)
ADD_SUBDIRECTORY(src-thyra)
#ADD_SUBDIRECTORY(src-petsc)
ADD_SUBDIRECTORY(src-loca)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(examples)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
