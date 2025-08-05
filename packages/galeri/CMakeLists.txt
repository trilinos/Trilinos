

#
# A) Define the package
#

TRIBITS_PACKAGE(Galeri)

#
# B) Set up package-specific options
#

TRIBITS_ADD_DEBUG_OPTION()

#
# C) Add the libraries, tests, and examples
#

#this dir only had "EXTRA_DIST" files in it
#ADD_SUBDIRECTORY(src-fem)

IF (${PACKAGE_NAME}_ENABLE_Xpetra)
  ADD_SUBDIRECTORY(src-xpetra)
ENDIF()


TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example-xpetra example-tpetra)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_FILES(
  src-pfem
  doc/AdvDiffSquare.png
  doc/L.*.png
  example-fem/TwoSquares.cpp
  src-fem/Galeri_FileGrid.h
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
