
# Declare all the packages that we use in tests here.  All these
# packages must be listed in either TEST_REQUIRED_DEP_PACKAGES,
# TEST_OPTIONAL_DEP_PACKAGES, LIB_REQUIRED_DEP_PACKAGES, or
# LIB_OPTIONAL_DEP_PACKAGES, in anasazi/cmake/Dependencies.cmake.  If
# any package in an ASSERT_DEFINED statement below is not listed in
# that file, the ASSERT_DEFINED will fail during package
# configuration.
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Belos)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Galeri)


if (${PACKAGE_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_EXECUTABLE(
    Epetra_BKS_Galeri_example
    SOURCES BlockKrylovSchurEpetra.cpp
    COMM serial mpi
  )
ENDIF ()


TRIBITS_ADD_EXECUTABLE(
  Epetra_BKS_example
  SOURCES BlockKrylovSchurEpetraEx.cpp
  COMM serial mpi
)


TRIBITS_ADD_EXECUTABLE(
  Epetra_BKS_SVD_example
  SOURCES BlockKrylovSchurEpetraExSVD.cpp
  COMM serial mpi
)






