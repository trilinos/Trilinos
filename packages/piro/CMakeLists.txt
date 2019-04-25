
#
# A) Define the package
#

TRIBITS_PACKAGE(Piro DISABLE_STRONG_WARNINGS)

#
# B) Set up package-specific options
#

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters)
GLOBAL_SET(HAVE_PIRO_THYRAEPETRAADAPTERS ${${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters)
GLOBAL_SET(HAVE_PIRO_THYRAEPETRAEXTADAPTERS ${${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Epetra)
GLOBAL_SET(HAVE_PIRO_EPETRA ${${PACKAGE_NAME}_ENABLE_Epetra})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
GLOBAL_SET(HAVE_PIRO_EPETRAEXT ${${PACKAGE_NAME}_ENABLE_EpetraExt})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)
GLOBAL_SET(HAVE_PIRO_TPETRA ${${PACKAGE_NAME}_ENABLE_Tpetra})

ASSERT_DEFINED(NOX_ENABLE_Epetra)
IF (${PACKAGE_NAME}_ENABLE_Epetra AND NOT NOX_ENABLE_Epetra)
  MESSAGE(FATAL_ERROR "Piro configuration error: Enabling Epetra in Piro also requires that Epetra be enabled in NOX (-DNOX_ENABLE_Epetra:BOOL=ON)")
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_OptiPack)
  MESSAGE(WARNING "***** WARNING:  OptiPack is deprecated in April 2019; it will be  removed in May 2019 ***** ")
  MESSAGE(WARNING "***** WARNING:  Option Piro_ENABLE_OptiPack is deprecated in April 2019; it will be removed in May 2019 ***** ")
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

# All tests are based on Epetra, currently
IF ( Piro_ENABLE_Epetra)
  TRIBITS_ADD_TEST_DIRECTORIES(test)
ENDIF ()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
