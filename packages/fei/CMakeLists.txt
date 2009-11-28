
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(FEI)

#
# B) Set up package-specific options
#

IF (Trilinos_ENABLE_Epetra)
ADD_OPTION_AND_DEFINE(
  FEI_EPETRA
  HAVE_FEI_EPETRA
  "FEI support for Epetra"
  ${Trilinos_ENABLE_Epetra}
)
ENDIF()

IF (Trilinos_ENABLE_AztecOO)
ADD_OPTION_AND_DEFINE(
  FEI_AZTECOO
  HAVE_FEI_AZTECOO
  "FEI support for AztecOO"
  ${Trilinos_ENABLE_AztecOO}
)
ENDIF()

IF (Trilinos_ENABLE_Amesos)
  ADD_OPTION_AND_DEFINE(
    FEI_AMESOS
    HAVE_FEI_AMESOS
    "FEI support for Amesos"
    ${Trilinos_ENABLE_Amesos}
  )
ENDIF()

IF (Trilinos_ENABLE_Tpetra)
ADD_OPTION_AND_DEFINE(
  FEI_TPETRA
  HAVE_FEI_TPETRA
  "FEI support for Tpetra"
  ${Trilinos_ENABLE_Tpetra}
)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(base)
ADD_SUBDIRECTORY(support-Trilinos)
ADD_SUBDIRECTORY(test_utils)


PACKAGE_ADD_TEST_DIRECTORIES(test)
PACKAGE_ADD_TEST_DIRECTORIES(unit_tests_base)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
