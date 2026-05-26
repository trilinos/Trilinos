

#
# A) Define the package
#

TRIBITS_PACKAGE(Galeri)

#
# B) Set up package-specific options
#

# Do explicit template instantiation (ETI) and testing over the same
# Scalar, LocalOrdinal, and GlobalOrdinal types as for Tpetra objects.
# However, by default, exclude all of Tpetra's GlobalOrdinal (GO)
# types from the ETI list.  This is because Tpetra adds all GO ETI
# types to the list of Scalar ETI types.  (This lets it implement
# certain communication patterns involving e.g., Vector<GO, LO, GO>.)
#

# Make sure that Tpetra actually defined these variables, even if they
# are empty.
ASSERT_DEFINED(TpetraCore_ETI_SCALARS_NO_ORDS)
ASSERT_DEFINED(TpetraCore_ETI_LORDS)
ASSERT_DEFINED(TpetraCore_ETI_GORDS)
ASSERT_DEFINED(TpetraCore_ETI_NODES)

# Promote these variables to be visible outside of the "scope" of this
# directory -- e.g., to other packages.
GLOBAL_SET(Galeri_ETI_SCALARS ${TpetraCore_ETI_SCALARS_NO_ORDS})
GLOBAL_SET(Galeri_ETI_LORDS   ${TpetraCore_ETI_LORDS})
GLOBAL_SET(Galeri_ETI_GORDS   ${TpetraCore_ETI_GORDS})
GLOBAL_SET(Galeri_ETI_NODES   ${TpetraCore_ETI_NODES})


TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

# we use this for testing as well as library eti, so do it regardless
TRIBITS_ADD_ETI_SUPPORT()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example-xpetra example-tpetra)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_FILES(
  src-pfem
  example-pfem
  doc/AdvDiffSquare.png
  doc/L.*.png
  example-fem/TwoSquares.cpp
  src-fem/Galeri_FileGrid.h
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
