

#
# A) Define the package
#

TRIBITS_PACKAGE(FEI)

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# B) Set up package-specific options
#

ASSERT_DEFINED(Trilinos_ENABLE_Epetra)
IF (NOT ${Trilinos_ENABLE_Epetra})
  MESSAGE(FATAL_ERROR "FEI requires Epetra be enabled")
ENDIF()
SET(FEI_EPETRA ON CACHE BOOL "FEI support for Epetra")
GLOBAL_SET(HAVE_FEI_EPETRA ON)

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
