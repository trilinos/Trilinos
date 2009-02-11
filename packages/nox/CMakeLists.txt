INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

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

IF(NOX_ENABLE_Sacado)
  SET(HAVE_LOCA_ANASAZI ON)
ENDIF()

###################################################
# Optional TPL Support
###################################################

###################################################
# Concrete implementations of NOX::Abstract Classes
###################################################

# LAPACK
SET(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK OFF CACHE BOOL "Enables support for LAPACK implementation of Abstract Group and Vector classes.  Defaults to OFF")

# Epetra - requires Epetra, AztecOO, and Ifpack
IF(NOX_ENABLE_AztecOO AND NOX_ENABLE_Epetra AND NOX_ENABLE_Ifpack)
  SET(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA_DEFAULT ON)
ELSE()
  SET(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA_DEFAULT OFF)
ENDIF()

SET(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA 
  ${NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA_DEFAULT} 
  CACHE BOOL "Enables support for Epetra implementation of Abstract Group and Vector classes.  Defaults to ON if Epetra, AztecOO, and Ifpack are enabled, otherwise defaults to OFF.")

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA)
  
  IF (NOT NOX_ENABLE_AztecOO OR NOT NOX_ENABLE_Epetra OR 
      NOT NOX_ENABLE_Ifpack)
    MESSAGE(FATAL_ERROR "ERROR: The NOX-Epetra Support has been enabled by the flag NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA but this requires the Trilinos libraries AztecOO, Epetra, and Ifpack to be enabled!  Please fix your configure script to enable all of these packages.  Current status is:\n Epetra:  ${NOX_ENABLE_Epetra} \n AztecOO: ${NOX_ENABLE_AztecOO} \n Ifpack:  ${NOX_ENABLE_Ifpack}")
  ENDIF()

ENDIF()

# Thyra
SET(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA ${NOX_ENABLE_Thyra} CACHE BOOL "Enables support for Thyra implementation of Abstract Group and Vector classes.  Defaults to ON if Thyra is enabled, otherwise defaults to OFF")

###################################################
###################################################

#
# C) Add the libraries, tests, and examples
#

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
