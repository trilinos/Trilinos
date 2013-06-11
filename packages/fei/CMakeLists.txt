
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(FEI)

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# B) Set up package-specific options
#

IF (${PACKAGE_NAME}_ENABLE_Epetra)
  TRIBITS_ADD_OPTION_AND_DEFINE(
    FEI_EPETRA
    HAVE_FEI_EPETRA
    "FEI support for Epetra"
    ${${PACKAGE_NAME}_ENABLE_Epetra}
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_AztecOO)
  TRIBITS_ADD_OPTION_AND_DEFINE(
    FEI_AZTECOO
    HAVE_FEI_AZTECOO
    "FEI support for AztecOO"
    ${${PACKAGE_NAME}_ENABLE_AztecOO}
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Belos)
  TRIBITS_ADD_OPTION_AND_DEFINE(
    FEI_BELOS
    HAVE_FEI_BELOS
    "FEI support for Belos"
    ${${PACKAGE_NAME}_ENABLE_Belos}
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Amesos)
    TRIBITS_ADD_OPTION_AND_DEFINE(
      FEI_AMESOS
      HAVE_FEI_AMESOS
      "FEI support for Amesos"
      ${${PACKAGE_NAME}_ENABLE_Amesos}
    )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_OPTION_AND_DEFINE(
    FEI_TPETRA
    HAVE_FEI_TPETRA
    "FEI support for Tpetra"
    ${${PACKAGE_NAME}_ENABLE_Tpetra}
  )
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(base)
ADD_SUBDIRECTORY(support-Trilinos)
ADD_SUBDIRECTORY(test_utils)


TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_TEST_DIRECTORIES(unit_tests_base)
TRIBITS_ADD_TEST_DIRECTORIES(unit_tests)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
