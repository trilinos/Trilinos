
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE( Anasazi ENABLE_SHADOWING_WARNINGS )

# 2009/05/01: rabartl: I have removed 'CLEANED' from above from now because it
# fails to compile with gcc 3.4.6 (see bug 4395).

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

PACKAGE_ADD_TEST_DIRECTORIES(test)

IF (Anasazi_ENABLE_Thyra)
  ADD_SUBDIRECTORY(thyra)
ENDIF()


IF (Anasazi_ENABLE_Tpetra)
  ADD_SUBDIRECTORY(tpetra)
ENDIF()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
