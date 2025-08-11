include(CombinedOption)

#
# A) Define the package
#

TRIBITS_PACKAGE(Piro DISABLE_STRONG_WARNINGS)

#
# B) Set up package-specific options
#





ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)
GLOBAL_SET(HAVE_PIRO_TPETRA ${${PACKAGE_NAME}_ENABLE_Tpetra})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Teko)
GLOBAL_SET(HAVE_PIRO_TEKO ${${PACKAGE_NAME}_ENABLE_Teko})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Stratimikos)
GLOBAL_SET(HAVE_PIRO_STRATIMIKOS ${${PACKAGE_NAME}_ENABLE_Stratimikos})


COMBINED_OPTION( PIRO_HAVE_EPETRA_STACK
  DEP_OPTIONS_NAMES ${PACKAGE_NAME}_ENABLE_Epetra
                    ${PACKAGE_NAME}_ENABLE_EpetraExt
                    ${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters
                    ${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters
                    ${PACKAGE_NAME}_ENABLE_Stratimikos
  DOCSTR "Enables support for Epetra objects. Defaults to ON if Epetra, EpetraExt, "
         "ThyraEpetraAdapters, ThyraEpetraExtAdapters and Stratimikos are enabled, "
         "otherwise defaults to OFF."
  )
MESSAGE(STATUS "PIRO_HAVE_EPETRA_STACK: ${PIRO_HAVE_EPETRA_STACK}")

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
