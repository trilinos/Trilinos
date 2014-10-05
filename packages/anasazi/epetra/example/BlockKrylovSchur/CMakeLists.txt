INCLUDE(TribitsAddExecutableAndTest)

# Declare all the packages that we use in tests here.  All these
# packages must be listed in either TEST_REQUIRED_DEP_PACKAGES,
# TEST_OPTIONAL_DEP_PACKAGES, LIB_REQUIRED_DEP_PACKAGES, or
# LIB_OPTIONAL_DEP_PACKAGES, in anasazi/cmake/Dependencies.cmake.  If
# any package in an ASSERT_DEFINED statement below is not listed in
# that file, the ASSERT_DEFINED will fail during package
# configuration.
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_AztecOO)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Belos)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Galeri)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)


if (${PACKAGE_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_EXECUTABLE(
    BlockKrylovSchurEpetra
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


IF (${PACKAGE_NAME}_ENABLE_EpetraExt)
  IF (${PACKAGE_NAME}_ENABLE_Triutils)
    TRIBITS_ADD_EXECUTABLE(
      Epetra_BKS_File_example
      SOURCES BlockKrylovSchurEpetraExFile.cpp
      COMM serial mpi
      )
  ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
ENDIF(${PACKAGE_NAME}_ENABLE_EpetraExt)


IF (${PACKAGE_NAME}_ENABLE_Amesos)
  TRIBITS_ADD_EXECUTABLE(
    Epetra_BKS_Amesos_example
    SOURCES BlockKrylovSchurEpetraExGenAmesos.cpp
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_ENABLE_EpetraExt)
    TRIBITS_ADD_EXECUTABLE(
      Epetra_BKS_Buckling_example
      SOURCES BlockKrylovSchurEpetraExBuckling.cpp
      COMM serial mpi
      )
  ENDIF(${PACKAGE_NAME}_ENABLE_EpetraExt)
ENDIF(${PACKAGE_NAME}_ENABLE_Amesos)


IF (${PACKAGE_NAME}_ENABLE_Ifpack)
  IF (${PACKAGE_NAME}_ENABLE_AztecOO)
    TRIBITS_ADD_EXECUTABLE(
      Epetra_BKS_AztecOO_example
      SOURCES BlockKrylovSchurEpetraExGenAztecOO.cpp
      COMM serial mpi
      )
  ENDIF(${PACKAGE_NAME}_ENABLE_AztecOO)

  IF (${PACKAGE_NAME}_ENABLE_Belos)
    TRIBITS_ADD_EXECUTABLE(
      Epetra_BKS_Belos_example
      SOURCES BlockKrylovSchurEpetraExGenBelos.cpp
      COMM serial mpi
      )
  ENDIF(${PACKAGE_NAME}_ENABLE_Belos)
ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)
