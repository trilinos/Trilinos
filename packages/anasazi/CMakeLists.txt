
INCLUDE(TrilinosPackageHelpers)
INCLUDE(Trilinos_Add_Option)

#
# A) Define the package
#

TRILINOS_PACKAGE(Anasazi)

#
# B) Set up package-specific options
#

# ...

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)


IF (Anasazi_ENABLE_Epetra)
  ADD_SUBDIRECTORY(epetra)
ENDIF()

TRILINOS_PACKAGE_ADD_TEST_DIRECTORIES(test)

#TRILINOS_PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

IF (Anasazi_ENABLE_Thyra)
  ADD_SUBDIRECTORY(thyra)
ENDIF()


IF (Anasazi_ENABLE_Tpetra)
  ADD_SUBDIRECTORY(tpetra)
ENDIF()

#
# D) Do standard postprocessing
#

TRILINOS_PACKAGE_POSTPROCESS()
