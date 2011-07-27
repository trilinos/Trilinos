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

# Prerelease code
OPTION(${PACKAGE_NAME}_BUILD_PRERELEASE
  "Build experimental code"
  OFF )
IF(NOX_BUILD_PRERELEASE)
  SET(WITH_PRERELEASE TRUE)
ENDIF()

# LOCA support
OPTION(${PACKAGE_NAME}_ENABLE_LOCA
  "Enable parameter LOCA support"
  ON )

# MF support
SET(NOX_ENABLE_MF_default OFF)
IF(TPL_ENABLE_MF)
  SET(NOX_ENABLE_MF_default ON)
ENDIF()
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_MF
  HAVE_LOCA_MF
  "Enable Multifario (MF) support in LOCA"
  ${NOX_ENABLE_MF_default} )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
  NOX_TEUCHOS_TIME_MONITOR
  "Build with Teuchos TimeMonitors enabled."
  "${Trilinos_ENABLE_TEUCHOS_TIME_MONITOR}" )

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
  DEP_OPTIONS_NAMES NOX_ENABLE_ThyraCore 
  DOCSTR "Enables support for Thyra implementation of Abstract Group"
    " and Vector classes.  Defaults to ON if Thyra is enabled,"
    " otherwise defaults to OFF"
  )

COMBINED_OPTION( NOX_ENABLE_THYRA_EPETRA_ADAPTERS
  DEP_OPTIONS_NAMES NOX_ENABLE_ThyraCore NOX_ENABLE_ThyraEpetraAdapters NOX_ENABLE_ThyraEpetraExtAdapters
  DOCSTR "Enables support for the Thyra Epetra(Ext) adapters "
    " Defaults to ON if the Thyra Epetra(Ext) adapters are both enabled,"
    " otherwise defaults to OFF"
  )

# PETSc
COMBINED_OPTION( NOX_ENABLE_ABSTRACT_IMPLEMENTATION_PETSC
  DEP_OPTIONS_NAMES NOX_ENABLE_PETSC
  DOCSTR "Enables support for PETSc implementation of Abstract Group and"
    " Vector classes.  Defaults to ON if NOX_ENABLE_PETSC is"
    " enabled, otherwise defaults to OFF."
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
ADD_SUBDIRECTORY(src-petsc)
ADD_SUBDIRECTORY(src-loca)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(examples)

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

PACKAGE_ARCH_EXCLUDE_FILES(
  src-loca/python
  test/lapack/LOCA_python
  )

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
