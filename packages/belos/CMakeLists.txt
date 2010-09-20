INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Belos)

#
# B) Set up package-specific options
#

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

IF (Belos_ENABLE_Epetra)
  ADD_SUBDIRECTORY(epetra)
ENDIF()

IF (Belos_ENABLE_Tpetra)
  ADD_SUBDIRECTORY(tpetra)
ENDIF()

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Tpetra_Timers
  HAVE_BELOS_TPETRA_TIMERS
  "Add timers to Belos::MultiVecTraits for Tpetra multivectors."
  NO
  )


#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
